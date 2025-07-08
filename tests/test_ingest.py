"""Integration tests for Ingest Pipeline; isolated tests for observable output

These tests verify that various matrix file types can be extracted and
transformed, as expected by code that loads transformed data into MongoDB.

Test doubles are used for test speed and isolation.

PREREQUISITE
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
Also, before 'EXAMPLES' commands, run `cd tests`.

EXAMPLES

# Run all tests
pytest

# Run all tests and see print() output
pytest -s

# Run only tests in test_ingest.py
pytest test_ingest.py

# Run only the test "test_good_metadata_file" in this test module
pytest test_ingest.py -k 'test_good_metadata_file' -s

# Run all tests, using multiple CPUs
# Details: https://pypi.org/project/pytest-xdist/
pytest -n auto -s

# Run all tests, show code coverage metrics
pytest --cov=../ingest/

"""

import unittest
from unittest.mock import patch, MagicMock
from test_dense import mock_load_r_files
import os
import glob

from pymongo.errors import AutoReconnect, BulkWriteError
from test_expression_files import mock_expression_load
from mock_gcp import mock_storage_client, mock_storage_blob

import config
from ingest_pipeline import (
    create_parser,
    validate_arguments,
    IngestPipeline,
    exit_pipeline,
    run_ingest,
    get_action_from_args,
)
from expression_files.expression_files import GeneExpression


def mock_load(self, *args, **kwargs):
    """Enables overwriting normal function with this placeholder.
    Returning the arguments enables tests to verify that the code invokes
    this method with expected argument values.

    TODO:
    Integrate MongoDB emulator for faster, higher-coverage tests (SCP-2000)

    This will enable us to also verify (and thus cover) loading-code *outputs*,
    unlike here where we merely give a way to verify loading-code *inputs*.
    Doing so via integration tests will isolate us from implementation changes.
    """
    self.load_args = args[0]
    self.load_kwargs = kwargs


# Mock method that writes to database
IngestPipeline.load = mock_load
GeneExpression.load = mock_expression_load

# Mixpanel logging needs config object instantiated
# tests should bypass mongo check for input BSON values
os.environ["BYPASS_MONGO_WRITES"] = "yes"
# Initialize global variables]
config.init("5d276a50421aa9117c982845", "5dd5ae25421aa910a723a337")
# restore environment variable to unset state
del os.environ['BYPASS_MONGO_WRITES']


class IngestTestCase(unittest.TestCase):
    @staticmethod
    @patch("google.cloud.storage.Blob", side_effect=mock_storage_blob)
    @patch("google.cloud.storage.Client", side_effect=mock_storage_client)
    def execute_ingest(args, mock_storage_client, mock_storage_blob):

        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)
        if "differential_expression" in arguments:
            # DE may use metadata or cluster file for annots BUT
            # IngestPipeline initialization assumes a "cell_metadata_file"
            arguments["cell_metadata_file"] = arguments["annotation_file"]
            # IngestPipeline initialization expects "name" and not "cluster_name"
            arguments["name"] = arguments["cluster_name"]


        ingest = IngestPipeline(**arguments)

        status, status_cell_metadata = run_ingest(ingest, arguments, parsed_args)

        return ingest, arguments, status, status_cell_metadata

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_dense_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract, transform, and load dense matrices"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "gs://fake-bucket/tests/data/dense_matrix_19_genes_1000_cells.txt",
            "--matrix-file-type",
            "dense",
        ]
        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_local_dense_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract and transform local dense matrices"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            "--matrix-file-type",
            "dense",
        ]
        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_local_compressed_dense_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract and transform local dense matrices
        from compressed file in the same manner as uncompressed file
        """

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/dense_matrix_19_genes_100k_cells.txt.gz",
            "--matrix-file-type",
            "dense",
        ]
        self.execute_ingest(args)

    def test_empty_dense_file(self):
        """Ingest Pipeline should fail gracefully when an empty file is given"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/empty_file.txt",
            "--matrix-file-type",
            "dense",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_empty_mtx_file(self):
        """Ingest Pipeline should fail gracefully when an empty file is given"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/empty_file.txt",
            "--barcode-file",
            "../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_mtx_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract and transform MTX matrix bundles"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            "--barcode-file",
            "../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        ]
        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_unsorted_mtx_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract and transform unsorted MTX matrix bundles"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/AB_toy_data_toy.unsorted_mtx.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/AB_toy_data_toy.genes.tsv",
            "--barcode-file",
            "../tests/data/AB_toy_data_toy.barcodes.tsv",
        ]
        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_ingest_zipped_mtx_matrix(self, mock_check_unique_cells):
        """Ingest Pipeline should extract and transform MTX matrix bundles"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/mtx/AB_toy_data_toy.unsorted_mtx.mtx.gz",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            "--barcode-file",
            "../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        ]
        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_remote_mtx_bundles(self, mock_check_unique_cells):
        """Ingest Pipeline should handle MTX matrix files fetched from bucket"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "gs://fake-bucket/tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "gs://fake-bucket/tests/data/mtx/AB_toy_data_toy.genes.tsv",
            "--barcode-file",
            "gs://fake-bucket/tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        ]

        self.execute_ingest(args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_mtx_bundle_argument_validation(self, mock_check_unique_cells):
        """Omitting --gene-file and --barcode-file in MTX ingest should error"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "../tests/data/mtx/matrix.mtx",
            "--matrix-file-type",
            "mtx",
        ]

        self.assertRaises(ValueError, self.execute_ingest, args)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_r_files,
    )
    def test_r_file_dense(self, mock_check_unique_cells, mock_load):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--matrix-file",
            "../tests/data/r_format_text.txt",
            "--matrix-file-type",
            "dense",
        ]
        self.execute_ingest(args)

    def test_good_metadata_file(self):
        """Ingest Pipeline should succeed for properly formatted metadata file"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cell_metadata",
            "--cell-metadata-file",
            "../tests/data/metadata_example.txt",
            "--study-accession",
            "SCP123",
            "--ingest-cell-metadata",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        # TODO:
        # After integrating MongoDB emulator (SCP-2000), refactor this test to
        # verify that status is "0", not [0, None] as done below.  This peculiar
        # expected value is tested because we intercept the load() function,
        # because we lack a MongoDB emulator.
        self.assertEqual(len(status), 2)
        self.assertEqual(status[1], None)

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=AutoReconnect,
    )
    def test_exponential_back_off_expression_file(
        self, mock_check_unique_cells, mock_load
    ):
        """Ingest Pipeline should not succeed if mongo cannot connect after 5
        reconnection tries.
        """
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_expression",
            "--taxon-name",
            "Homo sapiens",
            "--taxon-common-name",
            "human",
            "--ncbi-taxid",
            "9606",
            "--genome-assembly-accession",
            "GCA_000001405.15",
            "--genome-annotation",
            "Ensembl 94",
            "--matrix-file",
            "gs://fake-bucket/tests/data/AB_toy_data_toy.matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "gs://fake-bucket/tests/data/AB_toy_data_toy.genes.tsv",
            "--barcode-file",
            "gs://fake-bucket/tests/data/AB_toy_data_toy.barcodes.tsv",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_bad_metadata_file(self):
        """Ingest Pipeline should not succeed for misformatted metadata file"""

        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cell_metadata",
            "--cell-metadata-file",
            "../tests/data/metadata_bad.txt",
            "--study-accession",
            "SCP123",
            "--ingest-cell-metadata",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_bad_metadata_file_contains_coordinates(self):
        """Ingest Pipeline should not succeed for metadata file containing
        coordinates
        """
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cell_metadata",
            "--cell-metadata-file",
            "../tests/data/metadata_has_coordinate_header.txt",
            "--study-accession",
            "SCP123",
            "--ingest-cell-metadata",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_good_cluster_file(self):
        """Ingest Pipeline should succeed for properly formatted cluster file"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cluster",
            "--cluster-file",
            "../tests/data/cluster_example.txt",
            "--ingest-cluster",
            "--name",
            "cluster1",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        # TODO:
        # After integrating MongoDB emulator (SCP-2000), refactor this test to
        # verify that status is "0", not [None] as done below.  This peculiar
        # expected value is tested because we intercept the load() function,
        # because we lack a MongoDB emulator.
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], None)

    def test_bad_cluster_file(self):
        """Ingest Pipeline should fail for misformatted cluster file"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cluster",
            "--cluster-file",
            "../tests/data/cluster_bad.txt",
            "--ingest-cluster",
            "--name",
            "cluster1",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_bad_cluster_missing_coordinate_file(self):
        """Ingest Pipeline should fail for missing coordinate in cluster file"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_cluster",
            "--cluster-file",
            "../tests/data/cluster_bad_missing_coordinate.txt",
            "--ingest-cluster",
            "--name",
            "cluster1",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    @patch("ingest_pipeline.IngestPipeline.load_subsample", return_value=0)
    def test_subsample(self, mock_load_subsample):
        """When cell values in cluster are present in cell metadata file ingest should succeed."""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_subsample",
            "--cluster-file",
            "../tests/data/cluster_example.txt",
            "--name",
            "custer1",
            "--cell-metadata-file",
            "../tests/data/metadata_example.txt",
            "--subsample",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)

    def test_get_cluster_query(self):
        """When subsampling AnnData files cluster name should be appended to query"""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_subsample",
            "--cluster-file",
            "../tests/data/good_subsample_cluster.csv",
            "--name",
            "custer1",
            "--cell-metadata-file",
            "../tests/data/test_cell_metadata.csv",
            "--subsample",
        ]
        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)
        ingest = IngestPipeline(**arguments)

        mock_metrics = MagicMock()
        mock_metrics.get_properties.return_value = {"fileType": "AnnData"}
        with patch("config.get_metric_properties", return_value=mock_metrics):
            query = ingest.get_cluster_query()
            expected_keys = ['study_id', 'study_file_id', 'name']
            self.assertEqual(expected_keys, list(query.keys()))

        mock_metrics.get_properties.return_value = {"fileType": "Cluster"}
        with patch("config.get_metric_properties", return_value=mock_metrics):
            query = ingest.get_cluster_query()
            expected_keys = ['study_id', 'study_file_id']
            self.assertEqual(expected_keys, list(query.keys()))

    @patch("ingest_pipeline.IngestPipeline.load_subsample", return_value=0)
    def test_subsample_no_cell_intersection(self, mock_load_subsample):
        """When cell values in cluster are not present in cell metadata file ingest should fail."""
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_subsample",
            "--cluster-file",
            "../tests/data/good_subsample_cluster.csv",
            "--name",
            "custer1",
            "--cell-metadata-file",
            "../tests/data/test_cell_metadata.csv",
            "--subsample",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        with self.assertRaises(SystemExit) as cm, self.assertRaises(ValueError):
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_extract_cluster_file_from_anndata(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_anndata",
            "--ingest-anndata",
            "--extract",
            "['cluster']",
            "--anndata-file",
            "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad",
            "--obsm-keys",
            "['X_tsne']",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)
        filename = 'h5ad_frag.cluster.X_tsne.tsv.gz'
        self.assertTrue(os.path.isfile(filename))

        # clean up cluster file
        try:
            os.remove(filename)
        except:
            print(f"Error while deleting file : {filename}")

    def test_invalid_obsm_key_for_anndata(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_anndata",
            "--ingest-anndata",
            "--extract",
            "['cluster']",
            "--anndata-file",
            "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad",
            "--obsm-keys",
            "['foo']",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        self.assertEqual(status[0], 1)

    def test_extract_metadata_file_from_anndata(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_anndata",
            "--ingest-anndata",
            "--extract",
            "['metadata']",
            "--anndata-file",
            "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)
        filename = 'h5ad_frag.metadata.tsv.gz'
        self.assertTrue(os.path.isfile(filename))

        # clean up metadata file
        try:
            os.remove(filename)
        except:
            print(f"Error while deleting file : {filename}")

    def test_extract_processed_matrix_from_anndata(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_anndata",
            "--ingest-anndata",
            "--extract",
            "['processed_expression']",
            "--anndata-file",
            "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad",
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)
        mtx_filename = 'h5ad_frag.matrix.processed.mtx.gz'
        self.assertTrue(os.path.isfile(mtx_filename))
        barcodes_filename = 'h5ad_frag.barcodes.processed.tsv.gz'
        self.assertTrue(os.path.isfile(barcodes_filename))
        features_filename = 'h5ad_frag.features.processed.tsv.gz'
        self.assertTrue(os.path.isfile(features_filename))

        # clean up mtx files
        files = [mtx_filename, barcodes_filename, features_filename]

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_execute_de_dense(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "differential_expression",
            "--annotation-file",
            "../tests/data/differential_expression/de_dense_metadata.tsv",
            "--annotation-name",
            "cell_type__ontology_label",
            "--annotation-type",
            "group",
            "--annotation-scope",
            "study",
            "--de-type",
            "rest",
            "--cluster-file",
            "../tests/data/differential_expression/de_dense_cluster.tsv",
            "--cluster-name",
            "dense_de_integration",
            "--matrix-file-path",
            "../tests/data/differential_expression/de_dense_matrix.tsv",
            "--matrix-file-type",
            "dense",
            "--study-accession",
            "SCP123",
            "--differential-expression"
        ]

        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)

        expected_file = "dense_de_integration--cell_type__ontology_label--cholinergic_neuron--study--wilcoxon.tsv"
        expected_output_match = (
            "dense_de_integration--cell_type__ontology_label--*--study--wilcoxon.tsv"
        )

        files = glob.glob(expected_output_match)

        self.assertIn(
            expected_file, files, "Expected filename not in found files list"
        )

        # clean up DE outputs
        output_wildcard_match = f"../tests/dense_de_integration--cell_type__ontology_label*.tsv"
        files = glob.glob(output_wildcard_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_execute_de_sparse(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "differential_expression",
            "--annotation-file",
            "../tests/data/differential_expression/sparse/sparsemini_metadata.txt",
            "--annotation-name",
            "cell_type__ontology_label",
            "--annotation-type",
            "group",
            "--annotation-scope",
            "study",
            "--de-type",
            "rest",
            "--cluster-file",
            "../tests/data/differential_expression/sparse/sparsemini_cluster.txt",
            "--cluster-name",
            "sparse_de_integration",
            "--matrix-file-path",
            "../tests/data/differential_expression/sparse/sparsemini_matrix.mtx",
            "--gene-file",
            "../tests/data/differential_expression/sparse/sparsemini_dup_gene_name.tsv",
            "--barcode-file",
            "../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv",
            "--matrix-file-type",
            "mtx",
            "--study-accession",
            "SCP123",
            "--differential-expression"
        ]

        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)

        expected_file = "sparse_de_integration--cell_type__ontology_label--fibroblast--study--wilcoxon.tsv"
        expected_output_match = (
            "sparse_de_integration--cell_type__ontology_label--*--study--wilcoxon.tsv"
        )

        files = glob.glob(expected_output_match)

        self.assertIn(
            expected_file, files, "Expected filename not in found files list"
        )

        # clean up DE outputs
        output_wildcard_match = f"../tests/sparse_de_integration--cell_type__ontology_label*.tsv"
        files = glob.glob(output_wildcard_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_execute_de_anndata(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "differential_expression",
            "--annotation-file",
            "../tests/data/anndata/compliant_liver_h5ad_frag.metadata.tsv.gz",
            "--annotation-name",
            "cell_type__ontology_label",
            "--annotation-type",
            "group",
            "--annotation-scope",
            "study",
            "--de-type",
            "rest",
            "--cluster-file",
            "../tests/data/anndata/compliant_liver_h5ad_frag.cluster.X_umap.tsv.gz",
            "--cluster-name",
            "umap",
            "--matrix-file-path",
            "../tests/data/anndata/compliant_liver.h5ad",
            "--matrix-file-type",
            "h5ad",
            "--raw-location",
            ".raw",
            "--study-accession",
            "SCP123",
            "--differential-expression"
        ]

        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)

        expected_file = "umap--cell_type__ontology_label--plasma_cell--study--wilcoxon.tsv"
        expected_output_match = (
            "umap--cell_type__ontology_label--*--study--wilcoxon.tsv"
        )

        files = glob.glob(expected_output_match)

        self.assertIn(
            expected_file, files, "Expected filename not in found files list"
        )

        # clean up DE outputs
        output_wildcard_match = f"../tests/umap--cell_type__ontology_label*.tsv"
        files = glob.glob(output_wildcard_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_ingest_dot_plot(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_dot_plot_genes",
            "--cell-metadata-file",
            "../tests/data/differential_expression/sparse/sparsemini_metadata.txt",
            "--cluster-group-id",
            "686d8982a374d4e7fadc89a6",
            "--cluster-file",
            "../tests/data/differential_expression/sparse/sparsemini_cluster.txt",
            "--matrix-file-path",
            "../tests/data/differential_expression/sparse/sparsemini_matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/differential_expression/sparse/sparsemini_dup_gene_name.tsv",
            "--barcode-file",
            "../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv"
        ]
        ingest, arguments, status, status_cell_metadata = self.execute_ingest(args)

        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], 0)

    def test_get_action_from_args(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_subsample",
            "--cluster-file",
            "../tests/data/good_subsample_cluster.csv",
            "--name",
            "cluster1",
            "--cell-metadata-file",
            "../tests/data/test_cell_metadata.csv",
            "--subsample",
        ]
        self.assertEqual("ingest_subsample", get_action_from_args(args))
        bad_args = ["foo", "bar", "bing"]
        self.assertEqual("", get_action_from_args(bad_args))

    @patch("mongo_connection.MongoConnection.MAX_AUTO_RECONNECT_ATTEMPTS", 3)
    def test_insert_reconnect(self):
        args = [
            "--study-id",
            "5d276a50421aa9117c982845",
            "--study-file-id",
            "5dd5ae25421aa910a723a337",
            "ingest_subsample",
            "--cluster-file",
            "../tests/data/good_subsample_cluster.csv",
            "--name",
            "cluster1",
            "--cell-metadata-file",
            "../tests/data/test_cell_metadata.csv",
            "--subsample",
        ]
        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)
        ingest = IngestPipeline(**arguments)
        client_mock = MagicMock()
        ingest.db = client_mock
        docs = [
            {
                'id': 1,
                'name': 'foo',
                'study_id': 1,
                'study_file_id': 1,
                'array_index': 0,
                'linear_data_type': 'Cluster',
            },
            {
                'id': 2,
                'name': 'bar',
                'study_id': 1,
                'study_file_id': 1,
                'array_index': 0,
                'linear_data_type': 'Cluster',
            },
        ]
        ingest.insert_many("data_arrays", docs)
        client_mock["data_arrays"].insert_many.assert_called_with(docs)

        client_mock["data_arrays"].insert_many.side_effect = ValueError("Foo")
        self.assertRaises(Exception, ingest.insert_many, "data_arrays", docs)
        client_mock.reset_mock()

        # Test exponential back off for auto reconnect
        client_mock["data_arrays"].insert_many.side_effect = AutoReconnect
        self.assertRaises(AutoReconnect, ingest.insert_many, "data_arrays", docs)
        self.assertEqual(client_mock["data_arrays"].insert_many.call_count, 3)
        client_mock.reset_mock()

        def raiseError(*args, **kwargs):
            details = {
                "writeErrors": [
                    {
                        "code": 11000,
                        "op": {
                            'id': 1,
                            'name': 'foo',
                            'study_id': 1,
                            'study_file_id': 1,
                            'array_index': 0,
                            'linear_data_type': 'Cluster',
                        },
                    }
                ]
            }
            raise BulkWriteError(details)

        # Test exponential back off for BulkWriteError
        client_mock["data_arrays"].insert_many.side_effect = raiseError
        self.assertRaises(BulkWriteError, ingest.insert_many, "data_arrays", docs)
        self.assertEqual(client_mock["data_arrays"].insert_many.call_count, 3)


if __name__ == "__main__":
    unittest.main()
