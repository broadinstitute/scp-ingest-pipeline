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
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/valid_no_array_v2.0.0.txt --study-accession SCP123 --ingest-cell-metadata

# Ingest Cell Metadata file against convention
!! Please note that you must have a pre-configured BigQuery table available
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/valid_no_array_v2.0.0.txt --study-accession SCP123 --ingest-cell-metadata --validate-convention --bq-dataset cell_metadata --bq-table alexandria_convention

# Ingest dense file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --ncbi-taxid 9606 --matrix-file ../tests/data/dense_matrix_19_genes_1000_cells.txt --matrix-file-type dense

# Ingest loom file
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --matrix-file ../tests/data/test_loom.loom  --matrix-file-type loom --taxon-name 'Homo Sapiens' --taxon-common-name human

# Subsample cluster and metadata file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_subsample --cluster-file ../tests/data/test_1k_cluster_Data.csv --name custer1 --cell-metadata-file ../tests/data/test_1k_metadata_Data.csv --subsample

# Ingest mtx files
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --matrix-file ../tests/data/mtx/matrix.mtx --matrix-file-type mtx --gene-file ../tests/data/genes.tsv --barcode-file ../tests/data/barcodes.tsv
"""
import json
import logging
import os
import re
import sys
from contextlib import nullcontext
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401


from bson.objectid import ObjectId


# from google.cloud.logging.resource import Resource
try:
    # Used when importing internally and in tests
    from ingest_files import IngestFiles

    # For Mixpanel logging
    import config
    from monitoring.mixpanel_log import custom_metric
    from monitoring.metrics_service import MetricsService

    # For tracing
    from opencensus.ext.stackdriver.trace_exporter import StackdriverExporter
    from opencensus.trace.samplers import AlwaysOnSampler
    from opencensus.trace.tracer import Tracer
    from pymongo import MongoClient
    from subsample import SubSample
    from validation.validate_metadata import (
        report_issues,
        validate_input_metadata,
        write_metadata_to_bq,
    )
    from cell_metadata import CellMetadata
    from cli_parser import create_parser, validate_arguments
    from clusters import Clusters
    from expression_files.mtx import MTXIngestor
    from expression_files.dense_ingestor import DenseIngestor
    from monitor import setup_logger, log_exception

except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles
    from . import config
    from .monitoring.metrics_service import MetricsService
    from .subsample import SubSample
    from .monitoring.mixpanel_log import custom_metric
    from .validation.validate_metadata import (
        validate_input_metadata,
        report_issues,
        write_metadata_to_bq,
    )
    from .monitor import setup_logger, log_exception
    from .cell_metadata import CellMetadata
    from .clusters import Clusters
    from .expression_files.dense_ingestor import DenseIngestor
    from .expression_files.mtx import MTXIngestor
    from .cli_parser import create_parser, validate_arguments


class IngestPipeline:
    # File location for metadata json convention
    JSON_CONVENTION = (
        "../schema/alexandria_convention/alexandria_convention_schema.json"
    )

    # Logger provides more details for trouble shooting
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    user_logger = setup_logger(__name__ + ".usr", "user_log.txt", level=logging.ERROR)

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
        if os.environ.get("DATABASE_HOST") is not None:
            # Needed to run tests in CircleCI.
            # TODO (SCP-2000): Integrate MongoDB emulator to test_ingest.py, then remove this
            self.db = self.get_mongo_db()
        else:
            self.db = None
        self.cluster_file = cluster_file
        self.kwargs = kwargs
        self.cell_metadata_file = cell_metadata_file
        if "GOOGLE_CLOUD_PROJECT" in os.environ:
            # instantiate trace exporter
            exporter = StackdriverExporter(
                project_id=os.environ["GOOGLE_CLOUD_PROJECT"]
            )
            self.tracer = Tracer(exporter=exporter, sampler=AlwaysOnSampler())

        else:
            self.tracer = nullcontext()
        if ingest_cell_metadata:
            self.cell_metadata = self.initialize_file_connection(
                "cell_metadata", cell_metadata_file
            )
        if ingest_cluster:
            self.cluster = self.initialize_file_connection("cluster", cluster_file)
        if subsample:
            self.cluster_file = cluster_file
            self.cell_metadata_file = cell_metadata_file

    # Will be replaced by MongoConnection as defined in SCP-2629
    def get_mongo_db(self):
        host = os.environ["DATABASE_HOST"]
        user = os.environ["MONGODB_USERNAME"]
        password = os.environ["MONGODB_PASSWORD"]
        db_name = os.environ["DATABASE_NAME"]
        client = MongoClient(
            host,
            username=user,
            password=password,
            authSource=db_name,
            authMechanism="SCRAM-SHA-1",
        )

        return client[db_name]

    def initialize_file_connection(self, file_type, file_path):
        """Initializes connection to file.

            Returns:
                File object.
        """
        file_connections = {"cell_metadata": CellMetadata, "cluster": Clusters}

        return file_connections.get(file_type)(
            file_path,
            self.study_id,
            self.study_file_id,
            tracer=self.tracer,
            **self.kwargs,
        )

    def insert_many(self, collection_name, documents):
        self.db[collection_name].insert_many(documents)

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
                collection_name == "cell_metadata"
                and model["name"] == "NAME"
                and model["annotation_type"] == "TYPE"
            ):
                linear_id = ObjectId(self.study_id)
            else:
                linear_id = self.db[collection_name].insert_one(model).inserted_id
            for data_array_model in set_data_array_fn(
                linear_id, *set_data_array_fn_args, **set_data_array_fn_kwargs
            ):
                documents.append(data_array_model)
            # only insert documents if present
            if len(documents) > 0:
                self.insert_many("data_arrays", documents)
        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    def load_subsample(
        self, parent_collection_name, subsampled_data, set_data_array_fn, scope
    ):
        """Loads subsampled data into MongoDB"""
        documents = []
        try:
            for key_value in subsampled_data[0].items():
                annot_name = subsampled_data[1][0]
                annot_type = subsampled_data[1][1]
                sample_size = subsampled_data[2]
                query = {
                    "study_id": ObjectId(self.study_id),
                    "study_file_id": ObjectId(self.study_file_id),
                }
                # Query mongo for linear_id and 'name' of parent
                # Then return 'name' and 'id' fields from query results
                parent_data = self.db[parent_collection_name].find_one(
                    query, {"name": 1}
                )
                for model in set_data_array_fn(
                    (
                        key_value[0],  # NAMES, x, y, or z
                        # Cluster name provided from parent
                        parent_data["name"],
                        key_value[1],  # Subsampled data/values
                        ObjectId(self.study_file_id),
                        ObjectId(self.study_id),
                        parent_data["_id"],
                    ),
                    {
                        "subsample_annotation": f"{annot_name}--{annot_type}--{scope}",
                        "subsample_threshold": sample_size,
                    },
                ):
                    documents.append(model)
            self.db["data_arrays"].insert_many(documents)

        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    def upload_metadata_to_bq(self):
        """Uploads metadata to BigQuery"""
        if self.kwargs["validate_convention"] is not None:
            if (
                self.kwargs["validate_convention"]
                and self.kwargs["bq_dataset"]
                and self.kwargs["bq_table"]
            ):
                write_status = write_metadata_to_bq(
                    self.cell_metadata,
                    self.kwargs["bq_dataset"],
                    self.kwargs["bq_table"],
                )
                return write_status
            else:
                IngestPipeline.dev_logger.error(
                    "Erroneous call to upload_metadata_to_bq"
                )
                return 1
        return 0

    @custom_metric(config.get_metric_properties)
    def ingest_expression(self) -> int:
        """
        Ingests expression files.
        """
        self.expression_ingestor = None
        try:
            if MTXIngestor.matches_file_type(self.matrix_file_type):
                self.expression_ingestor = MTXIngestor(
                    self.matrix_file, self.study_id, self.study_file_id, **self.kwargs
                )
            if DenseIngestor.matches_file_type(self.matrix_file_type):
                self.expression_ingestor = DenseIngestor(
                    self.matrix_file,
                    self.study_id,
                    self.study_file_id,
                    tracer=self.tracer,
                    **self.kwargs,
                )
            self.expression_ingestor.execute_ingest()
        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    # More work needs to be done to fully remove ingest from IngestPipeline
    # Tracked in SCP-3023
    @custom_metric(config.get_metric_properties)
    def ingest_cell_metadata(self):
        """Ingests cell metadata files into Firestore."""
        validate_against_convention = False
        if self.kwargs["validate_convention"] is not None:
            if self.kwargs["validate_convention"]:
                validate_against_convention = True
        self.cell_metadata.preprocess(validate_against_convention)
        if self.cell_metadata.validate(validate_against_convention):
            IngestPipeline.dev_logger.info("Cell metadata file format valid")
            # Check file against metadata convention
            if validate_against_convention:
                if self.cell_metadata.conforms_to_metadata_convention():
                    IngestPipeline.dev_logger.info(
                        "Cell metadata file conforms to metadata convention"
                    )
                    pass
                else:
                    return 1

            for metadata_model in self.cell_metadata.execute_ingest():
                IngestPipeline.dev_logger.info(
                    f"Attempting to load cell metadata header : {metadata_model.annot_header}"
                )
                status = self.load(
                    self.cell_metadata.COLLECTION_NAME,
                    metadata_model.model,
                    self.cell_metadata.set_data_array,
                    metadata_model.annot_header,
                )
                if status != 0:
                    IngestPipeline.user_logger.error(
                        f"Loading cell metadata header : {metadata_model.annot_header} failed. Exiting program"
                    )
                    return status
            return status if status is not None else 1
        else:
            report_issues(self.cell_metadata)
            IngestPipeline.user_logger.error("Cell metadata file format invalid")
            return 1

    @custom_metric(config.get_metric_properties)
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
            report_issues(self.cluster)
            IngestPipeline.user_logger.error("Cluster file format invalid")
            return 1
        return status

    @custom_metric(config.get_metric_properties)
    def subsample(self):

        """Method for subsampling cluster and metadata files"""
        subsample = SubSample(
            cluster_file=self.cluster_file, cell_metadata_file=self.cell_metadata_file
        )
        for data in subsample.subsample("cluster"):
            load_status = self.load_subsample(
                Clusters.COLLECTION_NAME, data, subsample.set_data_array, "cluster"
            )

            if load_status != 0:
                return load_status

        if self.cell_metadata_file is not None:
            try:
                subsample.prepare_cell_metadata()
                # Get cell names from cluster and metadata files
                # strip of whitespace that pandas might add
                cluster_cell_names = map(
                    lambda s: s.strip(), SubSample.get_cell_names(subsample.file)
                )
                metadata_cell_names = map(
                    lambda s: s.strip(),
                    SubSample.get_cell_names(subsample.cell_metadata.file),
                )
                # Check that cell names in cluster file exist in cell metadata file
                if SubSample.has_cells_in_metadata_file(
                    metadata_cell_names, cluster_cell_names
                ):
                    for data in subsample.subsample("study"):
                        load_status = self.load_subsample(
                            Clusters.COLLECTION_NAME,
                            data,
                            subsample.set_data_array,
                            "study",
                        )
                        if load_status != 0:
                            return load_status
                else:
                    raise ValueError(
                        "Cluster file has cell names that are not present in cell metadata file."
                    )
            except Exception as e:
                log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
                return 1
        return 0


def run_ingest(ingest, arguments, parsed_args):
    """Runs Ingest Pipeline as indicated by CLI or importing (test) module
    """
    status = []
    status_cell_metadata = None
    # TODO: Add validation for gene file types
    if "matrix_file" in arguments:
        config.set_parent_event_name("ingest-pipeline:expression:ingest")
        status.append(ingest.ingest_expression())
    elif "ingest_cell_metadata" in arguments:
        if arguments["ingest_cell_metadata"]:
            config.set_parent_event_name("ingest-pipeline:cell_metadata:ingest")
            status_cell_metadata = ingest.ingest_cell_metadata()
            status.append(status_cell_metadata)
            if parsed_args.bq_table is not None and status_cell_metadata == 0:
                status_metadata_bq = ingest.upload_metadata_to_bq()
                status.append(status_metadata_bq)
    elif "ingest_cluster" in arguments:
        if arguments["ingest_cluster"]:
            config.set_parent_event_name("ingest-pipeline:cluster:ingest")
            status.append(ingest.ingest_cluster())
    elif "subsample" in arguments:
        if arguments["subsample"]:
            config.set_parent_event_name("ingest-pipeline:subsample:ingest")
            status_subsample = ingest.subsample()
            status.append(status_subsample)

    return status, status_cell_metadata


def exit_pipeline(ingest, status, status_cell_metadata, arguments):
    """Logs any errors, then exits Ingest Pipeline with standard OS code
    """
    if len(status) > 0:
        if all(i < 1 for i in status):
            sys.exit(os.EX_OK)
        else:
            # delocalize errors file
            for argument in list(arguments.keys()):
                captured_argument = re.match("(\w*file)$", argument)
                if captured_argument is not None:
                    study_file_id = arguments["study_file_id"]
                    matched_argument = captured_argument.groups()[0]
                    file_path = arguments[matched_argument]
                    if IngestFiles.is_remote_file(file_path):
                        # Delocalize support log
                        IngestFiles.delocalize_file(
                            study_file_id,
                            arguments["study_id"],
                            file_path,
                            "log.txt",
                            f"parse_logs/{study_file_id}/log.txt",
                        )
                        # Delocalize user log
                        IngestFiles.delocalize_file(
                            study_file_id,
                            arguments["study_id"],
                            file_path,
                            "user_log.txt",
                            f"parse_logs/{study_file_id}/user_log.txt",
                        )
                    # Need 1 argument that has a path to identify google bucket
                    # Break after first argument
                    break
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
    parsed_args = create_parser().parse_args()
    validate_arguments(parsed_args)
    arguments = vars(parsed_args)
    # Initialize global variables for current ingest job
    config.init(
        arguments["study_id"],
        arguments["study_file_id"],
        arguments["user_metrics_uuid"],
    )
    ingest = IngestPipeline(**arguments)
    status, status_cell_metadata = run_ingest(ingest, arguments, parsed_args)
    # Log Mixpanel events
    MetricsService.log(config.get_parent_event_name(), config.get_metric_properties())
    # Exit pipeline
    exit_pipeline(ingest, status, status_cell_metadata, arguments)


if __name__ == "__main__":
    main()
