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
import random
import sys
import unittest
from unittest.mock import patch

from gcp_mocks import mock_storage_client, mock_storage_blob

from google.cloud import firestore
import google.auth.credentials

sys.path.append("../ingest")
from ingest_pipeline import create_parser, validate_arguments, IngestPipeline


def _make_credentials():
    return unittest.mock.Mock(spec=google.auth.credentials.Credentials)


def get_random_study_accession():
    '''Randomly seeds to avoid collisions in Firestore emulator instance'''
    study_number = random.randint(1, 1_000_001)
    return f'SCP{study_number}'  # E.g. SCP1324


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

        study_accession = get_random_study_accession()

        args = [
            "--study-accession",
            study_accession,
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
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        query_results = [doc.to_dict() for doc in stream]

        num_results = len(query_results)
        self.assertEqual(num_results, 19)

    def test_ingest_missing_file(self):
        """Ingest Pipeline should throw error for missing file
        """

        study_accession = get_random_study_accession()

        args = [
            "--study-accession",
            study_accession,
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
            "gs://fake-bucket/remote-matrix-file-does-not-exist.txt",
            "--matrix-file-type",
            "dense",
        ]

        self.assertRaises(OSError, self.setup_ingest, args)

    def test_ingest_local_dense_matrix(self):
        """Ingest Pipeline should extract and transform local dense matrices
        """

        study_accession = get_random_study_accession()

        args = [
            "--study-accession",
            study_accession,
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
            "../tests/data/dense_matrix_19_genes_100k_cells.txt",
            "--matrix-file-type",
            "dense",
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        query_results = [doc.to_dict() for doc in stream]

        num_results = len(query_results)
        self.assertEqual(num_results, 19)

    def test_ingest_missing_local_file(self):
        """Ingest Pipeline should throw error for missing local file
        """

        study_accession = get_random_study_accession()

        args = [
            "--study-accession",
            study_accession,
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
            "--matrix-file /this/file/does/not_exist.txt",
            "--matrix-file-type",
            "dense",
        ]

        self.assertRaises(OSError, self.setup_ingest, args)

    def test_ingest_mtx_matrix(self):
        """Ingest Pipeline should extract and transform MTX matrix bundles
        """

        study_accession = get_random_study_accession()

        args = [
            "--study-accession",
            study_accession,
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
            "../tests/data/matrix.mtx",
            "--matrix-file-type",
            "mtx",
            "--gene-file",
            "../tests/data/genes.tsv",
            "--barcode-file",
            "../tests/data/barcodes.tsv",
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        query_results = [doc.to_dict() for doc in stream]

        num_results = len(query_results)
        self.assertEqual(num_results, 25)

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

    def test_mtx_bundle_argument_validation(self):
        """Omitting --gene-file and --barcode-file in MTX ingest should error
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
            "../tests/data/matrix.mtx",
            "--matrix-file-type",
            "mtx",
        ]

        self.assertRaises(ValueError, self.setup_ingest, args)


if __name__ == "__main__":
    unittest.main()
