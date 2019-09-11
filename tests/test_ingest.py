"""Integration tests for Ingest Pipeline; isolated tests for observable output

These tests verify that various matrix file types can be extracted and
transformed, as expected by code that loads transformed data into Firestore.

Mocks are used for test speed and isolation.

PREREQUISITE
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
Also, before "EXAMPLES" commands, run `cd tests`.

EXAMPLES

# Run all tests
pytest

# Run tests with names containing the string "dense"
pytest -k "dense"

# Run all tests in a manner that shows any print() statements
python3 test_ingest.py

# Run all tests, using multiple CPUs
# Details: https://pypi.org/project/pytest-xdist/
pytest -n auto

# Run all tests, show code coverage metrics
# For explanation of coverage report, see:
# https://coverage.readthedocs.io/en/v4.5.x/branch.html
coverage run --branch test_ingest.py; coverage report -m --include *scp-ingest-pipeline/ingest*

"""
import ast
import sys
import unittest
from unittest.mock import patch

from gcp_mocks import mock_storage_client, mock_storage_blob

from google.cloud import firestore
import google.auth.credentials

sys.path.append("../ingest")
from ingest_pipeline import create_parser, validate_arguments, IngestPipeline


# def mock_load_expression_data(self, *args, **kwargs):
#     """Enables overwriting normal function with this placeholder.
#     Returning the arguments enables tests to verify that the code invokes
#     this method with expected argument values.

#     TODO:
#     Use Google's official Firestore mock,
#     https://github.com/googleapis/google-cloud-ruby/blob/master/google-cloud-firestore/EMULATOR.md#google-cloud-firestore-emulator

#     This will enable us to also verify (and thus cover) loading code *outputs*,
#     unlike here where we merely give a way to verify loading code *inputs*.
#     Doing so via integration tests will isolate us from implementation changes.
#     """
#     self.load_expression_data_args = args
#     self.load_expression_data_kwargs = kwargs


def _make_credentials():
    return unittest.mock.Mock(spec=google.auth.credentials.Credentials)


def get_nth_gene_models(n, models, mock_dir):
    """Return Nth actual and expected gene models, using actual and mock data
    """
    # TODO: Dense loads models as a `list`, Mtx loads models as a `dict_values`
    # It seems both would ideally load using the same type.  Reconcile.
    if isinstance(models, list):
        # For Dense
        actual_model = models[n].__dict__
    else:
        # For Mtx
        actual_model = list(models)[n].__dict__

    # Uncomment to print out new baseline data
    # Process to update baselines is manual: copy and paste it into new file
    # TODO: Automate when reasonable
    print("actual_model")
    print(actual_model)

    with open(f"mock_data/{mock_dir}/gene_model_{n}.txt") as f:
        # Create a dictionary from the string-literal mock
        expected_model = ast.literal_eval(f.read())

    return actual_model, expected_model


# Mock method that loads data to Firestore
# IngestPipeline.load_expression_data = mock_load_expression_data


class IngestTestCase(unittest.TestCase):
    @patch("google.cloud.storage.Blob", side_effect=mock_storage_blob)
    @patch("google.cloud.storage.Client", side_effect=mock_storage_client)
    def setup_ingest(self, args, mock_storage_client, mock_storage_blob):

        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        credentials = _make_credentials()
        db = firestore.Client(project='test-project', credentials=credentials)
        arguments['db'] = db
        print('arguments')
        print(arguments)

        ingest = IngestPipeline(**arguments)

        if "matrix_file" in arguments:
            ingest.ingest_expression()
        elif "cell_metadata_file" in arguments:
            ingest.ingest_cell_metadata()
        elif "cluster_file" in arguments:
            ingest.ingest_cluster()

        return ingest

    def test_ingest_dense_matrix(self):
        """Ingest Pipeline should extract and transform dense matrices
        """

        args = [
            "--study-accession",
            "SCP1",
            "--file-id",
            "1234abc",
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
            "gs://fake-bucket/tests/data/dense_matrix_19_genes_100k_cells.txt",
            "--matrix-file-type",
            "dense",
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = genes.where(u'taxon_common_name', u'==', u'human').stream()

        query_results = [doc.to_dict() for doc in stream]

        num_results = len(query_results)
        self.assertEqual(num_results, 19)

    # def test_ingest_missing_file(self):
    #     """Ingest Pipeline should throw error for missing file
    #     """

    #     args = [
    #         "--study-accession",
    #         "SCP1",
    #         "--file-id",
    #         "1234abc",
    #         "ingest_expression",
    #         "--taxon-name",
    #         "Homo sapiens",
    #         "--taxon-common-name",
    #         "human",
    #         "--ncbi-taxid",
    #         "9606",
    #         "--genome-assembly-accession",
    #         "GCA_000001405.15",
    #         "--genome-annotation",
    #         "Ensembl 94",
    #         "--matrix-file",
    #         "gs://fake-bucket/remote-matrix-file-does-not-exist.txt",
    #         "--matrix-file-type",
    #         "dense",
    #     ]

    #     self.assertRaises(OSError, self.setup_ingest, args)

    # def test_ingest_local_dense_matrix(self):
    #     """Ingest Pipeline should extract and transform local dense matrices
    #     """

    #     args = [
    #         "--study-accession",
    #         "SCP1",
    #         "--file-id",
    #         "1234abc",
    #         "ingest_expression",
    #         "--taxon-name",
    #         "Homo sapiens",
    #         "--taxon-common-name",
    #         "human",
    #         "--ncbi-taxid",
    #         "9606",
    #         "--genome-assembly-accession",
    #         "GCA_000001405.15",
    #         "--genome-annotation",
    #         "Ensembl 94",
    #         "--matrix-file",
    #         "../tests/data/dense_matrix_19_genes_100k_cells.txt",
    #         "--matrix-file-type",
    #         "dense",
    #     ]
    #     ingest = self.setup_ingest(args)

    #     models = ingest.load_expression_data_args[0]

    #     # Verify that 19 gene models were passed into load method
    #     num_models = len(models)
    #     expected_num_models = 19
    #     self.assertEqual(num_models, expected_num_models)

    #     # Verify that the first gene model looks as expected
    #     mock_dir = "dense_matrix_19_genes_100k_cells_txt"
    #     model, expected_model = get_nth_gene_models(0, models, mock_dir)

    #     self.assertEqual(model, expected_model)

    # def test_ingest_missing_local_file(self):
    #     """Ingest Pipeline should throw error for missing local file
    #     """

    #     args = [
    #         "--study-accession",
    #         "SCP1",
    #         "--file-id",
    #         "1234abc",
    #         "ingest_expression",
    #         "--taxon-name",
    #         "Homo sapiens",
    #         "--taxon-common-name",
    #         "human",
    #         "--ncbi-taxid",
    #         "9606",
    #         "--genome-assembly-accession",
    #         "GCA_000001405.15",
    #         "--genome-annotation",
    #         "Ensembl 94",
    #         "--matrix-file",
    #         "--matrix-file /this/file/does/not_exist.txt",
    #         "--matrix-file-type",
    #         "dense",
    #     ]

    #     self.assertRaises(OSError, self.setup_ingest, args)

    # def test_ingest_mtx_matrix(self):
    #     """Ingest Pipeline should extract and transform MTX matrix bundles
    #     """

    #     args = [
    #         "--study-accession",
    #         "SCP1",
    #         "--file-id",
    #         "1234abc",
    #         "ingest_expression",
    #         "--taxon-name",
    #         "Homo sapiens",
    #         "--taxon-common-name",
    #         "human",
    #         "--ncbi-taxid",
    #         "9606",
    #         "--genome-assembly-accession",
    #         "GCA_000001405.15",
    #         "--genome-annotation",
    #         "Ensembl 94",
    #         "--matrix-file",
    #         "../tests/data/matrix.mtx",
    #         "--matrix-file-type",
    #         "mtx",
    #         "--gene-file",
    #         "../tests/data/genes.tsv",
    #         "--barcode-file",
    #         "../tests/data/barcodes.tsv",
    #     ]
    #     ingest = self.setup_ingest(args)

    #     models = ingest.load_expression_data_args[0]

    #     # Verify that 25 gene models were passed into load method
    #     num_models = len(models)
    #     expected_num_models = 25
    #     self.assertEqual(num_models, expected_num_models)

    #     # Verify that the first gene model looks as expected
    #     mock_dir = "matrix_mtx"
    #     model, expected_model = get_nth_gene_models(0, models, mock_dir)
    #     print(f"\n\n\n{model}\n\n\n")
    #     print(f"\n\n\n{expected_model}\n\n\n")
    #     self.assertEqual(model, expected_model)

    # def test_mtx_bundle_argument_validation(self):
    #     """Omitting --gene-file and --barcode-file in MTX ingest should error
    #     """

    #     args = [
    #         "--study-accession",
    #         "SCP1",
    #         "--file-id",
    #         "1234abc",
    #         "ingest_expression",
    #         "--taxon-name",
    #         "Homo sapiens",
    #         "--taxon-common-name",
    #         "human",
    #         "--ncbi-taxid",
    #         "9606",
    #         "--genome-assembly-accession",
    #         "GCA_000001405.15",
    #         "--genome-annotation",
    #         "Ensembl 94",
    #         "--matrix-file",
    #         "../tests/data/matrix.mtx",
    #         "--matrix-file-type",
    #         "mtx",
    #     ]

    #     self.assertRaises(ValueError, self.setup_ingest, args)


if __name__ == "__main__":
    unittest.main()
