"""Ingest Pipeline for ingesting expression, metadata and cluster
files into MongoDB.

DESCRIPTION
This CLI extracts and transforms different file types then writes them into
a remote MongoDB instance.

PREREQUISITES
See https://github.com/broadinstitute/scp-ingest-pipeline#prerequisites

EXAMPLES
# Takes expression file and stores it into MongoDB

# Ingest cluster file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cluster --cluster-file ../tests/data/test_1k_cluster_Data.csv --ingest-cluster --name cluster1 --domain-ranges "{'x':[-1, 1], 'y':[-1, 1], 'z':[-1, 1]}"

# Ingest Cell Metadata file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/valid_no_array_v1.1.3.tsv --study-accession SCP123 --ingest-cell-metadata

# Ingest Cell Metadata file against convention
!! Please note that you must have a pre-configured BigQuery table available
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/valid_no_array_v1.1.3.tsv --study-accession SCP123 --ingest-cell-metadata --validate-convention --bq-dataset cell_metadata --bq-table alexandria_convention

# Ingest dense file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --ncbi-taxid 9606 --matrix-file ../tests/data/dense_matrix_19_genes_1000_cells.txt --matrix-file-type dense

# Ingest loom file
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --matrix-file ../tests/data/test_loom.loom  --matrix-file-type loom --taxon-name 'Homo Sapiens' --taxon-common-name humans

# Subsample cluster and metadata file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_subsample --cluster-file ../tests/data/test_1k_cluster_Data.csv --name custer1 --cell-metadata-file ../tests/data/test_1k_metadata_Data.csv --subsample

# Ingest mtx files
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name humans --matrix-file ../tests/data/matrix.mtx --matrix-file-type mtx --gene-file ../tests/data/genes.tsv --barcode-file ../tests/data/barcodes.tsv
"""
import datetime
import json
import logging
import os
import pprint
import re
import sys
from contextlib import nullcontext
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

# import google.cloud.logging
from bson.objectid import ObjectId
# For tracing
from opencensus.ext.stackdriver.trace_exporter import StackdriverExporter
from opencensus.trace.samplers import AlwaysOnSampler
from opencensus.trace.tracer import Tracer
from pymongo import InsertOne, MongoClient
from pymongo.errors import BulkWriteError

# from google.cloud.logging.resource import Resource
#
# try:
# Used when importing internally and in tests
from ingest_files import IngestFiles
from subsample import SubSample
from loom import Loom
from validation.validate_metadata import (
    validate_input_metadata,
    report_issues,
    write_metadata_to_bq,
)
from monitor import setup_logger, log, trace
from cell_metadata import CellMetadata
from clusters import Clusters
from expression_files.dense_command import DenseCommand
from expression_files.dense import Dense
from invoker import IngestInvoker
from mtx import Mtx
from cli_parser import create_parser, validate_arguments
# except ImportError:
#     # Used when importing as external package, e.g. imports in single_cell_portal code
#     from .ingest_files import IngestFiles
#     from .subsample import SubSample
#     from .loom import Loom
#     from .validation.validate_metadata import (
#         validate_input_metadata,
#         report_issues,
#         write_metadata_to_bq,
#     )
#     from .monitor import setup_logger, log, trace
#     from .cell_metadata import CellMetadata
#     from .clusters import Clusters
#     from .dense import Dense
#     from .mtx import Mtx
#     from .cli_parser import create_parser, validate_arguments


class IngestPipeline(object):
    # File location for metadata json convention
    JSON_CONVENTION = (
        '../schema/alexandria_convention/alexandria_convention_schema.json'
    )
    logger = logging.getLogger(__name__)
    error_logger = setup_logger(
        __name__ + '_errors', 'errors.txt', level=logging.ERROR)
    info_logger = setup_logger(__name__, 'info.txt')
    my_debug_logger = log(error_logger)

    def __init__(
        self,
        study_id: str,
        study_file_id: str,
        matrix_file: str = None,
        matrix_file_type: str = None,
        cell_metadata_file: str = None,
        cluster_file: str = None,
        subsample=False,
        ingest_cell_metadata=False,
        ingest_cluster=False,
        **kwargs,
    ):
        """Initializes variables in Ingest Pipeline"""
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.matrix_file = matrix_file
        self.matrix_file_type = matrix_file_type
        self.extra_log_params = {'study_id': self.study_id, 'duration': None}
        if os.environ.get('DATABASE_HOST') is not None:
            # Needed to run tests in CircleCI.  TODO: add mock, remove this
            self.error_logger.error(f"Database host is: {os.environ.get('DATABASE_HOST')}", extra=self.extra_log_params)
            print(f"Database host is: {os.environ.get('DATABASE_HOST')}")
            self.db = self.get_mongo_db()
        else:
            self.error_logger.error(f"Database host is: NONE", extra=self.extra_log_params)
            self.db = None
        self.cluster_file = cluster_file
        self.kwargs = kwargs
        self.cell_metadata_file = cell_metadata_file
        if 'GOOGLE_CLOUD_PROJECT' in os.environ:
            # instantiate trace exporter
            exporter = StackdriverExporter(
                project_id=os.environ['GOOGLE_CLOUD_PROJECT']
            )
            self.tracer = Tracer(exporter=exporter, sampler=AlwaysOnSampler())

        else:
            self.tracer = nullcontext()
        if ingest_cell_metadata:
            self.cell_metadata = self.initialize_file_connection(
                "cell_metadata", cell_metadata_file
            )
        if ingest_cluster:
            self.cluster = self.initialize_file_connection(
                "cluster", cluster_file)
        if matrix_file is None:
            self.matrix = matrix_file

        if subsample:
            self.cluster_file = cluster_file
            self.cell_metadata_file = cell_metadata_file

    @my_debug_logger()
    def get_mongo_db(self):
        host = os.environ['DATABASE_HOST']
        user = os.environ['MONGODB_USERNAME']
        password = os.environ['MONGODB_PASSWORD']
        db_name = os.environ['DATABASE_NAME']
        print(f'db_name:{db_name},  user:{user}, host:{host}')
        self.error_logger.error(f'db_name:{db_name},  user:{user}, host:{host}', extra=self.extra_log_params)
        client = MongoClient(
            host,
            username=user,
            password=password,
            authSource=db_name,
            authMechanism='SCRAM-SHA-1',
        )

        return client[db_name]

    def initialize_file_connection(self, file_type, file_path):
        """Initializes connection to file.

            Returns:
                File object.
        """
        self.error_logger.error(f'Trying to initialize file', extra=self.extra_log_params)
        print(f'Trying to initialize file')
        # Mtx file types not included because class declaration is different
        file_connections = {
            "dense": Dense,
            "cell_metadata": CellMetadata,
            "cluster": Clusters,
            "mtx": Mtx,
            "loom": Loom,
        }

        return file_connections.get(file_type)(
            file_path,
            self.study_id,
            self.study_file_id,
            tracer=self.tracer,
            **self.kwargs,
        )

    def close_matrix(self):
        """Closes connection to file"""
        self.matrix.close()

    # @profile
    # TODO: Make @profile conditional (SCP-2081)
    def insert_many(self, collection_name, documents):
        self.db[collection_name].insert_many(documents)

    def ingest_expression(self):
        expression_command = None
        print(self.matrix_file_type)
        if self.matrix_file_type == 'dense':
            dense = Dense(self.matrix_file,
                          self.study_id,
                          self.study_file_id,
                          tracer=self.tracer,
                          **self.kwargs,)  # Dense receiver
            expression_command = DenseCommand(dense)  # Dense concrete command
        invoker = IngestInvoker(expression_command)  # invoker
        invoker.execute()  # Starting the method calls

    @trace
    # @profile
    def load(
        self,
        collection_name,
        model,
        set_data_array_fn,
        *set_data_array_fn_args,
        **set_data_array_fn_kwargs,
    ):
        documents = []
        try:
            # hack to avoid inserting invalid CellMetadata object from first column
            # TODO: implement method similar to kwargs solution in ingest_expression
            if (
                collection_name == 'cell_metadata'
                and model['name'] == 'NAME'
                and model['annotation_type'] == 'TYPE'
            ):
                linear_id = ObjectId(self.study_id)
            else:
                linear_id = self.db[collection_name].insert_one(
                    model).inserted_id
            for data_array_model in set_data_array_fn(
                linear_id, *set_data_array_fn_args, **set_data_array_fn_kwargs
            ):
                documents.append(data_array_model)
            # only insert documents if present
            if len(documents) > 0:
                self.insert_many('data_arrays', documents)
        except Exception as e:
            self.error_logger.error(e, extra=self.extra_log_params)
            return 1
        return 0

    # @my_debug_logger()
    def ingest_cell_metadata(self):
        """Ingests cell metadata files into Firestore."""
        self.cell_metadata.preprocess()
        if self.cell_metadata.validate():
            self.info_logger.info(
                f'Cell metadata file format valid', extra=self.extra_log_params
            )
            # Check file against metadata convention
            if self.kwargs['validate_convention'] is not None:
                if self.kwargs['validate_convention']:
                    if self.conforms_to_metadata_convention():
                        self.info_logger.info(
                            f'Cell metadata file conforms to metadata convention',
                            extra=self.extra_log_params,
                        )
                        pass
                    else:
                        return 1

            self.cell_metadata.reset_file()
            for metadataModel in self.cell_metadata.transform():
                self.info_logger.info(
                    f'Attempting to load cell metadata header : {metadataModel.annot_header}',
                    extra=self.extra_log_params,
                )
                status = self.load(
                    self.cell_metadata.COLLECTION_NAME,
                    metadataModel.model,
                    self.cell_metadata.set_data_array,
                    metadataModel.annot_header,
                )
                if status != 0:
                    self.error_logger.error(
                        f'Loading cell metadata header : {metadataModel.annot_header} failed. Exiting program',
                        extra=self.extra_log_params,
                    )
                    return status
            return status if status is not None else 1
        else:
            report_issues(self.cell_metadata)
            self.error_logger.error(
                f'Cell metadata file format invalid', extra=self.extra_log_params
            )
            return 1

    @my_debug_logger()
    def ingest_cluster(self):
        """Ingests cluster files."""
        if self.cluster.validate():
            annotation_model = self.cluster.transform()
            status = self.load(
                self.cluster.COLLECTION_NAME,
                annotation_model,
                self.cluster.get_data_array_annot,
            )
            if status != 0:
                return status
        # Incorrect file format
        else:
            self.error_logger.error(
                f'Cluster file format invalid', extra=self.extra_log_params
            )
            return 1
        return status

    @my_debug_logger()
    def subsample(self):
        """Method for subsampling cluster and metadata files"""
        subsample = SubSample(
            cluster_file=self.cluster_file, cell_metadata_file=self.cell_metadata_file
        )
        for data in subsample.subsample('cluster'):
            load_status = self.load_subsample(
                Clusters.COLLECTION_NAME, data, subsample.set_data_array, 'cluster'
            )

            if load_status != 0:
                return load_status

        if self.cell_metadata_file is not None:
            try:
                subsample.prepare_cell_metadata()
                for data in subsample.subsample('study'):
                    load_status = self.load_subsample(
                        Clusters.COLLECTION_NAME,
                        data,
                        subsample.set_data_array,
                        'study',
                    )
                    if load_status != 0:
                        return load_status
            except Exception as e:
                self.error_logger.error(e, extra=self.extra_log_params)
                return 1
        return 0


def run_ingest(ingest, arguments, parsed_args):
    """Runs Ingest Pipeline as indicated by CLI or importing (test) module
    """
    status = []
    status_cell_metadata = None
    logging.basicConfig(filename='errors.txt', level=logging.DEBUG)
    logging.debug(f'In run_ingest: {arguments}')
    print(f'In run_ingest: {arguments}')
    # TODO: Add validation for gene file types
    if "matrix_file" in arguments:
        status.append(ingest.ingest_expression())
    elif "ingest_cell_metadata" in arguments:
        if arguments["ingest_cell_metadata"]:
            status_cell_metadata = ingest.ingest_cell_metadata()
            status.append(status_cell_metadata)
            if parsed_args.bq_table is not None and status_cell_metadata == 0:
                status_metadata_bq = ingest.upload_metadata_to_bq()
                status.append(status_metadata_bq)
    elif "ingest_cluster" in arguments:
        if arguments["ingest_cluster"]:
            status.append(ingest.ingest_cluster())
    elif "subsample" in arguments:
        if arguments["subsample"]:
            status_subsample = ingest.subsample()
            status.append(status_subsample)

    return status, status_cell_metadata


def exit_pipeline(ingest, status, status_cell_metadata, arguments):
    """Logs any errors, then exits Ingest Pipeline with standard OS code
    """
    logging.basicConfig(filename='errors.txt', level=logging.DEBUG)
    logging.debug("exit_pipeline now")
    print('exit_pipeline now')
    for argument in list(arguments.keys()):

        captured_argument = re.match("(\w*file)$", argument)
        if captured_argument is not None:
            study_file_id = arguments['study_file_id']
            matched_argument = captured_argument.groups()[0]
            file_path = arguments[matched_argument]
            if IngestFiles.is_remote_file(file_path):
                IngestFiles.delocalize_file(
                    study_file_id,
                    arguments['study_id'],
                    file_path,
                    'errors.txt',
                    f'parse_logs/{study_file_id}/errors.txt',
                )
            # Need 1 argument that has a path to identify google bucket
            # Break after first argument
            break
    if len(status) > 0:
        if all(i < 1 for i in status):
            sys.exit(os.EX_OK)
        else:
            # delocalize errors file
            # for argument in list(arguments.keys()):
            #     captured_argument = re.match("(\w*file)$", argument)
            #     if captured_argument is not None:
            #         study_file_id = arguments['study_file_id']
            #         matched_argument = captured_argument.groups()[0]
            #         file_path = arguments[matched_argument]
            #         if IngestFiles.is_remote_file(file_path):
            #             IngestFiles.delocalize_file(
            #                 study_file_id,
            #                 arguments['study_id'],
            #                 file_path,
            #                 'errors.txt',
            #                 f'parse_logs/{study_file_id}/errors.txt',
            #             )
                    # Need 1 argument that has a path to identify google bucket
                    # Break after first argument
                    # break
            if status_cell_metadata is not None:
                if status_cell_metadata > 0 and ingest.cell_metadata.is_remote_file:
                    # PAPI jobs failing metadata validation against convention report
                    # will have "unexpected exit status 65 was not ignored"
                    # EX_DATAERR (65) The input data was incorrect in some way.
                    # note that failure to load to MongoDB also triggers this error
                    sys.exit(os.EX_DATAERR)
            sys.exit(1)


def main() -> None:
    """Enables running Ingest Pipeline via CLI

    Args:
        None

    Returns:
        None
    """
    logging.basicConfig(filename='errors.txt', level=logging.DEBUG)
    print("Made it into main")
    logging.debug("Made it into main")
    parsed_args = create_parser().parse_args()
    logging.debug(f'Here are the parsed args: {parsed_args}')
    print(f'Here are the parsed args: {parsed_args}')
    validate_arguments(parsed_args)
    arguments = vars(parsed_args)
    logging.debug(f'Initailizing ingest pipeline')
    print('Initailizing ingest pipeline')
    ingest = IngestPipeline(**arguments)
    print('Starting ingest')
    logging.debug(f'Starting ingest')
    status, status_cell_metadata = run_ingest(ingest, arguments, parsed_args)
    exit_pipeline(ingest, status, status_cell_metadata, arguments)


if __name__ == "__main__":
    main()
