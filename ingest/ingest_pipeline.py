"""Ingest Pipeline for ingesting expression, metadata and cluster
files into MongoDB.

DESCRIPTION
This CLI extracts and transforms different file types then writes them into
a remote MongoDB instance.

PREREQUISITES
See https://github.com/broadinstitute/scp-ingest-pipeline#prerequisites

DEVELOPER SETUP (see README.md#Install and ../scripts/setup-mongo-dev.sh)

EXAMPLES
# Takes expression file and stores it into MongoDB

# Ingest cluster file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cluster --cluster-file ../tests/data/test_1k_cluster_data.csv --ingest-cluster --name cluster1 --domain-ranges "{'x':[-1, 1], 'y':[-1, 1], 'z':[-1, 1]}"

# Ingest Cell Metadata file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt --study-accession SCP123 --ingest-cell-metadata

# Ingest Cell Metadata file against convention
!! Please note that you must have a pre-configured BigQuery table available
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_cell_metadata --cell-metadata-file ../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt --study-accession SCP123 --ingest-cell-metadata --validate-convention

# Ingest Cell Metadata file against convention with has_<modality> metadata
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 ingest_cell_metadata --cell-metadata-file ../tests/data/annotation/metadata/convention/brain_rf1/patchseq_classic_metadata_has_modality_10.tsv --study-accession SCPPR344 --ingest-cell-metadata --validate-convention --has-modality "['electrophysiology', 'morphology']"

# Ingest dense file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --ncbi-taxid 9606 --matrix-file ../tests/data/dense_matrix_19_genes_1000_cells.txt --matrix-file-type dense

# Ingest AnnData file basic "does it open in Scanpy" validation
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/test.h5ad

# Ingest AnnData file - cluster extraction
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --extract "['cluster']" --ingest-anndata --anndata-file ../tests/data/anndata/test.h5ad --obsm-keys "['X_tsne']"

# Subsample cluster and metadata file
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_subsample --cluster-file ../tests/data/test_1k_cluster_data.csv --name cluster1 --cell-metadata-file ../tests/data/test_1k_metadata_Data.csv --subsample

# Ingest MTX files
python ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --matrix-file ../tests/data/mtx/matrix.mtx --matrix-file-type mtx --gene-file ../tests/data/genes.tsv --barcode-file ../tests/data/barcodes.tsv

# Ingest AnnData as reference file
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/test.h5ad

# Ingest AnnData - happy path cluster-only extraction
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad --extract "['cluster']" --obsm-keys "['X_umap','X_tsne']"

# Ingest AnnData - happy path metadata-only extraction
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad  --extract "['metadata']"

# Ingest AnnData - happy path processed expression data only extraction
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad  --extract "['processed_expression']"

# Ingest AnnData - happy path raw count cell name only extraction
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad  --extract "['raw_counts']" --raw-location ".raw"

# Ingest AnnData - happy path cluster and metadata extraction
python ingest_pipeline.py  --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 ingest_anndata --ingest-anndata --anndata-file ../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad  --extract "['cluster', 'metadata']" --obsm-keys "['X_umap','X_tsne']"

# Differential expression analysis (dense matrix)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression --annotation-name cell_type__ontology_label --annotation-type group --annotation-scope study --matrix-file-path ../tests/data/differential_expression/de_dense_matrix.tsv --matrix-file-type dense --annotation-file ../tests/data/differential_expression/de_dense_metadata.tsv --cluster-file ../tests/data/differential_expression/de_dense_cluster.tsv --cluster-name de_integration --study-accession SCPdev --differential-expression

# Differential expression analysis (sparse matrix)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression --annotation-name cell_type__ontology_label --annotation-type group --annotation-scope study --matrix-file-path ../tests/data/differential_expression/sparse/sparsemini_matrix.mtx --gene-file ../tests/data/differential_expression/sparse/sparsemini_features.tsv --barcode-file ../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv --matrix-file-type mtx --annotation-file ../tests/data/differential_expression/sparse/sparsemini_metadata.txt --cluster-file ../tests/data/differential_expression/sparse/sparsemini_cluster.txt --cluster-name de_sparse_integration --study-accession SCPsparsemini --differential-expression

# Differential expression analysis (h5ad matrix, raw count in raw slot)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression  --raw-location '.raw' --annotation-name cell_type__ontology_label --de-type rest  --annotation-type group --annotation-scope study --annotation-file ../tests/data/anndata/compliant_liver_h5ad_frag.metadata.tsv.gz --cluster-file ../tests/data/anndata/compliant_liver_h5ad_frag.cluster.X_umap.tsv.gz --cluster-name umap --matrix-file-path ../tests/data/anndata/compliant_liver.h5ad  --matrix-file-type h5ad --study-accession SCPdev --differential-expression

# Differential expression analysis (h5ad matrix, raw count in adata.layers['counts'])
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression  --raw-location 'counts' --annotation-name cell_type__ontology_label --de-type rest  --annotation-type group --annotation-scope study --annotation-file ../tests/data/anndata/compliant_liver_h5ad_frag.metadata.tsv.gz --cluster-file ../tests/data/anndata/compliant_liver_h5ad_frag.cluster.X_umap.tsv.gz --cluster-name umap --matrix-file-path ../tests/data/anndata/compliant_liver_layers_counts.h5ad  --matrix-file-type h5ad --study-accession SCPdev --differential-expression

# Pairwise differential expression analysis (dense matrix)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression --annotation-name cell_type__ontology_label --de-type pairwise --group1 "['cholinergic neuron']" --group2 "cranial somatomotor neuron" --annotation-type group --annotation-scope study --matrix-file-path ../tests/data/differential_expression/de_dense_matrix.tsv --matrix-file-type dense --annotation-file ../tests/data/differential_expression/de_dense_metadata.tsv --cluster-file ../tests/data/differential_expression/de_dense_cluster.tsv --cluster-name de_integration --study-accession SCPdev --differential-expression

# Pairwise differential expression analysis (sparse matrix)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression --annotation-name cell_type__ontology_label --de-type pairwise --group1 "['endothelial cell']" --group2 "smooth muscle cell" --annotation-type group --annotation-scope study --matrix-file-path ../tests/data/differential_expression/sparse/sparsemini_matrix.mtx --gene-file ../tests/data/differential_expression/sparse/sparsemini_features.tsv --barcode-file ../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv --matrix-file-type mtx --annotation-file ../tests/data/differential_expression/sparse/sparsemini_metadata.txt --cluster-file ../tests/data/differential_expression/sparse/sparsemini_cluster.txt --cluster-name de_sparse_integration --study-accession SCPsparsemini --differential-expression

# Pairwise differential expression analysis (h5ad matrix, raw count in raw slot)
python ingest_pipeline.py --study-id addedfeed000000000000000 --study-file-id dec0dedfeed1111111111111 differential_expression  --raw-location '.raw' --annotation-name cell_type__ontology_label --de-type pairwise --group1 "mature B cell" --group2 "plasma cell" --annotation-type group --annotation-scope study --annotation-file ../tests/data/anndata/compliant_liver_h5ad_frag.metadata.tsv.gz --cluster-file ../tests/data/anndata/compliant_liver_h5ad_frag.cluster.X_umap.tsv.gz --cluster-name umap --matrix-file-path ../tests/data/anndata/compliant_liver.h5ad  --matrix-file-type h5ad --study-accession SCPdev --differential-expression

"""

import json
import logging
import os
import re
import sys
import re
from contextlib import nullcontext
import traceback
from typing import Dict, Generator, List, Tuple, Union
from wsgiref.simple_server import WSGIRequestHandler  # noqa: F401
from bson.objectid import ObjectId


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
from mongo_connection import MongoConnection, graceful_auto_reconnect
from subsample import SubSample
from validation.validate_metadata import (
    report_issues,
    validate_input_metadata,
)
from cell_metadata import CellMetadata
from cli_parser import create_parser, validate_arguments
from clusters import Clusters
from expression_files.mtx import MTXIngestor
from expression_files.dense_ingestor import DenseIngestor
from monitor import setup_logger, log_exception
from de import DifferentialExpression
from author_de import AuthorDifferentialExpression
from expression_writer import ExpressionWriter
from rank_genes import RankGenes

# scanpy uses anndata python package, disamibguate local anndata
# using underscore https://peps.python.org/pep-0008/#naming-conventions
from anndata_ import AnnDataIngestor


class IngestPipeline:
    # File location for metadata json convention
    JSON_CONVENTION = (
        "../schema/alexandria_convention/alexandria_convention_schema.json"
    )

    # array of actions to use when reporting to Mixpanel
    ACTION_NAMES = [
        'ingest_cluster',
        'ingest_cell_metadata',
        'ingest_expression',
        'ingest_anndata',
        'ingest_subsample',
        'ingest_differential_expression',
        'differential_expression',
        'render_expression_arrays',
        'rank_genes',
    ]

    # Logger provides more details for trouble shooting
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")
    user_logger = setup_logger(__name__ + ".usr", "user_log.txt", level=logging.ERROR)

    def __init__(
        self,
        study_id: str,
        study_file_id: str,
        matrix_file: str = None,
        matrix_file_path: str = None,
        matrix_file_type: str = None,
        cell_metadata_file: str = None,
        cluster_file: str = None,
        anndata_file: str = None,
        subsample=False,
        validate_cell_metadata=False,
        ingest_cell_metadata=False,
        ingest_cluster=False,
        differential_expression=False,
        ingest_differential_expression=False,
        **kwargs,
    ):
        """Initializes variables in Ingest Pipeline"""
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.matrix_file = matrix_file
        self.matrix_file_path = matrix_file_path
        self.matrix_file_type = matrix_file_type
        if os.environ.get("DATABASE_HOST") is not None:
            # Needed to run tests in CircleCI.
            # TODO (SCP-2000): Integrate MongoDB emulator to test_ingest.py, then remove this
            self.db = self.get_mongo_db()
        else:
            self.db = None
        self.cluster_file = cluster_file
        self.anndata_file = anndata_file
        self.kwargs = kwargs
        self.cell_metadata_file = cell_metadata_file
        self.props = {}
        if "GOOGLE_CLOUD_PROJECT" in os.environ:
            # instantiate trace exporter
            exporter = StackdriverExporter(
                project_id=os.environ["GOOGLE_CLOUD_PROJECT"]
            )
            self.tracer = Tracer(exporter=exporter, sampler=AlwaysOnSampler())

        else:
            self.tracer = nullcontext()
        if ingest_cell_metadata or differential_expression:
            self.cell_metadata = self.initialize_file_connection(
                "cell_metadata", cell_metadata_file
            )
        if ingest_cluster or differential_expression:
            self.cluster = self.initialize_file_connection("cluster", cluster_file)
        if subsample:
            self.cluster_file = cluster_file
            self.cell_metadata_file = cell_metadata_file
        if ingest_differential_expression:
            self.cell_metadata = ""
            self.cluster = ""

    # Will be replaced by MongoConnection as defined in SCP-2629
    def get_mongo_db(self):
        return MongoConnection()._client

    def initialize_file_connection(self, file_type, file_path):
        """Initializes connection to file.

        Returns:
            File object.
        """
        file_connections = {"cell_metadata": CellMetadata, "cluster": Clusters}
        try:
            return file_connections.get(file_type)(
                file_path,
                self.study_id,
                self.study_file_id,
                tracer=self.tracer,
                **self.kwargs,
            )
        except ValueError as v:
            # Caution: recording errorTypes in this manner can clobber other collected errors.
            # ValueErrors during file connection indicate file cannot be processed
            # this logging approach should not lose collected file validation information
            if str(v).startswith("could not convert"):
                config.get_metric_properties().update(
                    {"errorTypes": ["content:type:not-numeric"]}
                )
            elif str(v).startswith("Unable to parse"):
                config.get_metric_properties().update(
                    {"errorTypes": ["format:cap:unique"]}
                )
            else:
                config.get_metric_properties().update(
                    {"errorTypes": ["parse:unhandled"]}
                )
            self.report_validation("failure")
            raise ValueError(v)

    @graceful_auto_reconnect
    def insert_many(self, collection_name, documents):
        if not config.bypass_mongo_writes():
            self.db[collection_name].insert_many(documents)

    def insert_one(self, collection_name, model):
        if not config.bypass_mongo_writes():
            linear_id = self.db[collection_name].insert_one(model).inserted_id
            return linear_id

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
                linear_id = self.insert_one(collection_name, model)
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
                query = self.get_cluster_query()

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
            self.insert_many("data_arrays", documents)

        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    def get_cluster_query(self):
        """Generate MongoDB query to load ClusterGroup to set association IDs when subsampling"""
        query = {
            "study_id": ObjectId(self.study_id),
            "study_file_id": ObjectId(self.study_file_id),
        }

        # if this is an AnnData file, we need to append in the cluster name, otherwise AnnData studies with
        # multiple clusters will fail subsampling as the first cluster is always returned from the query
        file_type = config.get_metric_properties().get_properties().get('fileType')
        if file_type and file_type == "AnnData":
            query["name"] = self.kwargs.get("name")

        return query

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
            self.report_validation("failure")
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        self.report_validation("success")
        return 0

    # More work needs to be done to fully remove ingest from IngestPipeline
    # Tracked in SCP-3023
    @custom_metric(config.get_metric_properties)
    def validate_cell_metadata(self):
        """If indicated, validate cell metadata against convention
        otherwise, just preprocess file
        """
        validate_against_convention = False
        if self.kwargs["validate_convention"] is not None:
            if self.kwargs["validate_convention"]:
                validate_against_convention = True
        self.cell_metadata.preprocess(validate_against_convention)
        if self.cell_metadata.validate(validate_against_convention):
            IngestPipeline.dev_logger.info("Cell metadata file structure valid")
            # Check file against metadata convention
            if validate_against_convention:
                if self.kwargs['has_modality'] is not None:
                    self.cell_metadata.booleanize_modality_metadata()
                if self.cell_metadata.conforms_to_metadata_convention():
                    IngestPipeline.dev_logger.info(
                        "Cell metadata file conforms to metadata convention"
                    )
                else:
                    config.get_metric_properties().update(self.cell_metadata.props)
                    self.report_validation("failure")
                    return 1
            self.report_validation("success")
        else:
            report_issues(self.cell_metadata)
            config.get_metric_properties().update(self.cell_metadata.props)
            self.report_validation("failure")
            IngestPipeline.user_logger.error("Cell metadata file format invalid")
            return 1
        return 0

    @custom_metric(config.get_metric_properties)
    def ingest_cell_metadata(self):
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

    @custom_metric(config.get_metric_properties)
    def ingest_cluster(self):
        """Ingests cluster files."""
        if self.cluster.validate():
            self.report_validation("success")
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
            config.get_metric_properties().update(self.cluster.props)
            self.report_validation("failure")
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
                    # Caution: recording errorTypes in this manner can clobber other collected errors.
                    # In subsampling, known failure modes are ValueErrors which stop processing so
                    # this logging approach should not lose file validation information
                    config.get_metric_properties().update(
                        {"errorTypes": ["content:missing:values-across-files"]}
                    )
                    self.report_validation("failure")
                    raise ValueError(
                        "Cluster file has cell names that are not present in cell metadata file."
                    )
            except Exception as e:
                # ToDo ingest.props["errorType"] = "subsample:"
                log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
                return 1
        return 0

    @custom_metric(config.get_metric_properties)
    def extract_from_anndata(self):
        """Extract data subsets from anndata file per SCP filetypes."""
        self.anndata = AnnDataIngestor(
            self.anndata_file, self.study_id, self.study_file_id, **self.kwargs
        )
        if self.anndata.basic_validation():
            # Get metadata extraction parameters and perform extraction
            if self.kwargs.get("extract") and "metadata" in self.kwargs.get("extract"):
                metadata_filename = "h5ad_frag.metadata.tsv"
                # TODO (SCP-5104): perform check for successful extraction or report failure and exit
                AnnDataIngestor.generate_metadata_file(
                    self.anndata.adata, metadata_filename
                )
            # Get cluster extraction parameters and perform extraction
            if self.kwargs.get("extract") and "cluster" in self.kwargs.get("extract"):
                if not self.kwargs["obsm_keys"]:
                    self.kwargs["obsm_keys"] = ["X_tsne"]
                # TODO (SCP-5104): perform check for successful extraction or report failure and exit
                try:
                    for key in self.kwargs["obsm_keys"]:
                        AnnDataIngestor.generate_cluster_header(self.anndata.adata, key)
                        AnnDataIngestor.generate_cluster_type_declaration(
                            self.anndata.adata, key
                        )
                        AnnDataIngestor.generate_cluster_body(self.anndata.adata, key)
                except KeyError as e:
                    msg = f"Unable to extract cluster data from anndata file. Please check the provided obsm key, {e}."
                    self.report_validation("failure")
                    log_exception(
                        IngestPipeline.dev_logger, IngestPipeline.user_logger, msg
                    )
                    return 1
                except Exception as e:
                    log_exception(
                        IngestPipeline.dev_logger, IngestPipeline.user_logger, e
                    )
                    return 1
            # process matrix data
            ### TODO (SCP-5102, SCP-5103): how to associate "raw_count" cells to anndata file
            if self.kwargs.get("extract") and "processed_expression" in self.kwargs.get(
                "extract"
            ):
                self.anndata.generate_processed_matrix(self.anndata.adata)

            if self.kwargs.get('extract') and "raw_counts" in self.kwargs.get(
                'extract'
            ):
                if self.anndata.validate_raw_location():
                    self.anndata.ingest_raw_cells()
                else:
                    self.report_validation("failure")
                    return 1
            self.report_validation("success")
            return 0
        # scanpy unable to open AnnData file
        else:
            self.report_validation("failure")
            return 1

    def calculate_de(self):
        """Run differential expression analysis"""
        try:
            de = DifferentialExpression(
                cluster=self.cluster,
                cell_metadata=self.cell_metadata,
                matrix_file_path=self.matrix_file_path,
                matrix_file_type=self.matrix_file_type,
                **self.kwargs,
            )
            de.execute_de()
        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        # ToDo: surface failed DE for analytics (SCP-4206)
        return 0

    def calculate_author_de(self):
        """Run differential expression analysis"""
        try:
            kwargs = self.kwargs
            if "comparison_group" not in kwargs:
                kwargs["comparison_group"] = None

            # Map canonical DE headers to headers in actual DE file to
            header_refmap = {
                "gene": kwargs["gene_header"],
                "group": kwargs["group_header"],
                "comparison_group": kwargs["comparison_group_header"],
                "size": kwargs["size_metric"],
                "significance": kwargs["significance_metric"],
            }

            author_de = AuthorDifferentialExpression(
                kwargs["cluster_name"],
                kwargs["annotation_name"],
                kwargs["study_accession"],
                kwargs["annotation_scope"],
                kwargs["method"],
                kwargs["differential_expression_file"],
                header_refmap,
            )
            author_de.execute()

            # execute function as MAIN then create method to search for the files created, delocalize them to bucket as with de.py. adjust other places it needs to be called
        except Exception as e:
            print(traceback.format_exc())
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        # ToDo: surface failed DE for analytics (SCP-4206)
        return 0

    def render_expression_arrays(self):
        try:
            exp_writer = ExpressionWriter(
                matrix_file_path=self.matrix_file_path,
                matrix_file_type=self.matrix_file_type,
                cluster_file_path=self.cluster_file,
                **self.kwargs,
            )
            exp_writer.render_artifacts()
        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    def rank_genes(self):
        try:
            kwargs = self.kwargs
            RankGenes(kwargs["study_accession"], kwargs["publication"])
        except Exception as e:
            log_exception(IngestPipeline.dev_logger, IngestPipeline.user_logger, e)
            return 1
        return 0

    def report_validation(self, status):
        self.props["status"] = status
        config.get_metric_properties().update(self.props)
        MetricsService.log("file-validation", config.get_metric_properties())


def run_ingest(ingest, arguments, parsed_args):
    """Runs Ingest Pipeline as indicated by CLI or importing (test) module"""
    status = []
    status_cell_metadata = None
    # TODO: Add validation for gene file types
    if "matrix_file" in arguments:
        config.set_parent_event_name("ingest-pipeline:expression:ingest")
        status.append(ingest.ingest_expression())
    elif "ingest_cell_metadata" in arguments:
        if arguments["ingest_cell_metadata"]:
            config.set_parent_event_name("ingest-pipeline:cell_metadata:ingest")
            status_cell_metadata_validation = ingest.validate_cell_metadata()
            status.append(status_cell_metadata_validation)
            if status_cell_metadata_validation == 0:
                if ingest.kwargs['has_modality'] is not None:
                    ingest.cell_metadata.file = CellMetadata.restore_modality_metadata(
                        ingest.cell_metadata
                    )
                status_cell_metadata_ingest = ingest.ingest_cell_metadata()
                status.append(status_cell_metadata_ingest)
    elif "ingest_cluster" in arguments:
        if arguments["ingest_cluster"]:
            config.set_parent_event_name("ingest-pipeline:cluster:ingest")
            status.append(ingest.ingest_cluster())
    elif "subsample" in arguments:
        if arguments["subsample"]:
            config.set_parent_event_name("ingest-pipeline:subsample:ingest")
            status_subsample = ingest.subsample()
            status.append(status_subsample)
    elif "ingest_anndata" in arguments:
        if arguments["ingest_anndata"]:
            config.set_parent_event_name("ingest-pipeline:anndata:ingest")
            status_anndata = ingest.extract_from_anndata()
            status.append(status_anndata)
    elif "differential_expression" in arguments:
        config.set_parent_event_name("ingest-pipeline:differential-expression")
        status_de = ingest.calculate_de()
        status.append(status_de)
        print(f"STATUS after DE {status}")
    elif "ingest_differential_expression" in arguments:
        config.set_parent_event_name("ingest-pipeline:differential-expression:ingest")
        status_de = ingest.calculate_author_de()
        status.append(status_de)
        print(f"STATUS after ingest author DE: {status}")

    elif "render_expression_arrays" in arguments:
        config.set_parent_event_name("image-pipeline:render-expression-arrays")
        status_exp_writer = ingest.render_expression_arrays()
        status.append(status_exp_writer)

    elif "rank_genes" in arguments:
        config.set_parent_event_name("image-pipeline:rank-genes")
        status_rank_genes = ingest.rank_genes()
        status.append(status_rank_genes)

    return status, status_cell_metadata


def get_delocalization_info(arguments):
    """extract info on study file for delocalization decision-making"""
    for argument in list(arguments.keys()):
        captured_argument = re.match("(\w*file)$", argument)
        if captured_argument is not None:
            study_file_id = arguments["study_file_id"]
            matched_argument = captured_argument.groups()[0]
            file_path = arguments[matched_argument]

            # Need 1 argument that has a path to identify google bucket
            # Break after first argument
            break
    return file_path, study_file_id


def exit_pipeline(ingest, status, status_cell_metadata, arguments):
    """Logs any errors, then exits Ingest Pipeline with standard OS code"""
    if len(status) > 0:
        # for successful DE jobs, need to delocalize results
        if (
            "differential_expression" in arguments
            or "ingest_differential_expression" in arguments
        ) and all(i < 1 for i in status):
            file_path, study_file_id = get_delocalization_info(arguments)
            print('file_path', file_path)
            print('study_file_id', study_file_id)
            # append status?
            if IngestFiles.is_remote_file(file_path):
                files_to_match = DifferentialExpression.string_for_output_match(
                    arguments
                )
                DifferentialExpression.delocalize_de_files(
                    file_path, study_file_id, files_to_match
                )
        # for successful anndata jobs, need to delocalize intermediate ingest files
        elif arguments.get("extract") and all(i < 1 for i in status):
            file_path, study_file_id = get_delocalization_info(arguments)
            # append status?
            files_to_delocalize = []
            if IngestFiles.is_remote_file(file_path):
                if "cluster" in arguments.get("extract"):
                    files_to_delocalize.extend(
                        AnnDataIngestor.clusterings_to_delocalize(arguments)
                    )
                if "metadata" in arguments.get("extract"):
                    metadata_filename = f"h5ad_frag.metadata.tsv.gz"
                    files_to_delocalize.append(metadata_filename)
                if "processed_expression" in arguments.get("extract"):
                    mtx = "h5ad_frag.matrix.processed.mtx.gz"
                    barcodes = "h5ad_frag.barcodes.processed.tsv.gz"
                    features = "h5ad_frag.features.processed.tsv.gz"
                    mtx_bundle = [mtx, barcodes, features]
                    files_to_delocalize.extend(mtx_bundle)
                AnnDataIngestor.delocalize_extracted_files(
                    file_path,
                    study_file_id,
                    arguments["study_accession"],
                    files_to_delocalize,
                )
        # all non-DE, non-anndata ingest jobs can exit on success
        elif all(i < 1 for i in status):
            sys.exit(os.EX_OK)
        else:
            if "rank_genes" in arguments:
                study_file_id = "rank_genes"
                file_path = f"gs://{arguments.get('bucket_name')}/blank"
            else:
                file_path, study_file_id = get_delocalization_info(arguments)
            if IngestFiles.is_remote_file(file_path):
                if "differential_expression" in arguments:
                    log_path = (
                        f"parse_logs/differential_expression/{study_file_id}/log.txt"
                    )
                else:
                    log_path = f"parse_logs/{study_file_id}/log.txt"
                # Delocalize support log
                IngestFiles.delocalize_file(file_path, "log.txt", log_path)
                # Delocalize user log
                user_log_path = f"parse_logs/{study_file_id}/user_log.txt"
                IngestFiles.delocalize_file(file_path, "user_log.txt", user_log_path)
            if status_cell_metadata is not None:
                if status_cell_metadata > 0 and ingest.cell_metadata.is_remote_file:
                    # PAPI jobs failing metadata validation against convention report
                    # will have "unexpected exit status 65 was not ignored"
                    # EX_DATAERR (65) The input data was incorrect in some way.
                    # note that failure to load to MongoDB also triggers this error
                    sys.exit(os.EX_DATAERR)
            sys.exit(1)


def get_action_from_args(arguments):
    """Get the action from list of arguments denoting which data type is being ingested/extracted"""
    action = list(set(IngestPipeline.ACTION_NAMES) & set(arguments))
    return action[0] if action else ""


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
    if "differential_expression" in arguments:
        # DE may use metadata or cluster file for annots BUT
        # IngestPipeline initialization assumes a "cell_metadata_file"
        arguments["cell_metadata_file"] = arguments["annotation_file"]
        # IngestPipeline initialization expects "name" and not "cluster_name"
        arguments["name"] = arguments["cluster_name"]

    # Initialize global variables for current ingest job
    config.init(
        arguments["study_id"],
        arguments["study_file_id"],
        arguments["user_metrics_uuid"],
        get_action_from_args(arguments),
    )
    ingest = IngestPipeline(**arguments)
    status, status_cell_metadata = run_ingest(ingest, arguments, parsed_args)
    # Print metrics properties
    metrics_dump = config.get_metric_properties().get_properties()
    for key in metrics_dump.keys():
        print(f"{key}: {metrics_dump[key]}")

    # Log Mixpanel events
    MetricsService.log(config.get_parent_event_name(), config.get_metric_properties())
    # Exit pipeline
    arguments["study_accession"] = metrics_dump["studyAccession"]
    exit_pipeline(ingest, status, status_cell_metadata, arguments)


if __name__ == "__main__":
    main()
