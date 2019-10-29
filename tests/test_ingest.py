"""Integration tests for Ingest Pipeline; isolated tests for observable output

These tests verify that various matrix file types can be extracted and
transformed, as loaded into (an official emulator of) Firestore

Test doubles are used for test speed and isolation.

PREREQUISITE
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
Also, before 'EXAMPLES' commands, run `cd tests`.

EXAMPLES

# Run all tests
pytest

# Run tests with names containing the string 'dense'
pytest -k 'dense'

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
import random
import sys
import unittest
from unittest.mock import patch

from gcp_mocks import mock_storage_client, mock_storage_blob

from google.cloud import firestore
import google.auth.credentials

sys.path.append('../ingest')
from ingest_pipeline import create_parser, validate_arguments, IngestPipeline


def _make_credentials():
    return unittest.mock.Mock(spec=google.auth.credentials.Credentials)


def get_random_study_accession():
    """Randomly seeds to avoid collisions in Firestore emulator instance"""
    study_number = random.randint(1, 1_000_001)
    return f'SCP{study_number}'  # E.g. SCP1324


def get_nth_gene_docs(n, docs, mock_dir):
    """Return Nth actual and expected gene documents, using actual and mock data
    """

    # Firestore does not return results in the same order every time,
    # so sort documents by name to enable comparing Nth query result.
    docs = sorted(docs, key=lambda doc: doc['name'])

    actual_doc = docs[n]

    # Uncomment to print out new baseline data
    # Process to update baselines is manual: copy and paste it into new file
    # TODO: Automate when reasonable
    # print(f'actual_doc: {actual_doc}')

    with open(f'mock_data/{mock_dir}/gene_doc_{n}.txt') as f:
        # Create a dictionary from the string-literal mock
        expected_doc = ast.literal_eval(f.read())

    # Expected study accession is different on every test run,
    # to avoid collisions in Firestore emulator instance.
    # So we remove this key from the actual and expected docs to avoid
    # false positive test failures.
    del actual_doc['study_accession']
    del expected_doc['study_accession']

    return actual_doc, expected_doc


class IngestTestCase(unittest.TestCase):
    @patch('google.cloud.storage.Blob', side_effect=mock_storage_blob)
    @patch('google.cloud.storage.Client', side_effect=mock_storage_client)
    def setup_ingest(self, args, mock_storage_client, mock_storage_blob):

        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        credentials = _make_credentials()
        db = firestore.Client(project='test-project', credentials=credentials)
        arguments['db'] = db

        ingest = IngestPipeline(**arguments)

        if 'matrix_file' in arguments:
            ingest.ingest_expression()
        elif 'ingest_cell_metadata' in arguments:
            if arguments['ingest_cell_metadata']:
                ingest.ingest_cell_metadata()
        elif 'ingest_cluster' in arguments:
            if arguments['ingest_cluster']:
                ingest.ingest_cluster()
        elif 'subsample' in arguments:
            if arguments['subsample']:
                ingest.subsample()
        return ingest

    def test_ingest_dense_matrix(self):
        """Ingest Pipeline should extract, transform, and load dense matrices
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            'gs://fake-bucket/tests/data/dense_matrix_19_genes_100k_cells.txt',
            '--matrix-file-type',
            'dense',
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        docs = [doc.to_dict() for doc in stream]

        # Verify that 19 gene docs were written to Firestore
        num_docs = len(docs)
        self.assertEqual(num_docs, 19)

        # Verify that the first gene document looks as expected
        mock_dir = 'dense_matrix_19_genes_100k_cells_txt'
        doc, expected_doc = get_nth_gene_docs(0, docs, mock_dir)
        self.assertEqual(doc, expected_doc)

    def test_ingest_missing_file(self):
        """Ingest Pipeline should throw error for missing file
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            'gs://fake-bucket/remote-matrix-file-does-not-exist.txt',
            '--matrix-file-type',
            'dense',
        ]

        self.assertRaises(OSError, self.setup_ingest, args)

    def test_ingest_local_dense_matrix(self):
        """Ingest Pipeline should extract, transform, and load local dense matrices
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/dense_matrix_19_genes_100k_cells.txt',
            '--matrix-file-type',
            'dense',
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        docs = [doc.to_dict() for doc in stream]

        num_docs = len(docs)
        self.assertEqual(num_docs, 19)

    def test_ingest_missing_local_file(self):
        """Ingest Pipeline should throw error for missing local file
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '--matrix-file /this/file/does/not_exist.txt',
            '--matrix-file-type',
            'dense',
        ]

        self.assertRaises(OSError, self.setup_ingest, args)

    def test_ingest_mtx_matrix(self):
        """Ingest Pipeline should extract, transform, and load MTX matrix bundles
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/matrix.mtx',
            '--matrix-file-type',
            'mtx',
            '--gene-file',
            '../tests/data/genes.tsv',
            '--barcode-file',
            '../tests/data/barcodes.tsv',
        ]
        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        docs = [doc.to_dict() for doc in stream]

        # Verify that 25 gene docs were written to Firestore
        num_docs = len(docs)
        self.assertEqual(num_docs, 25)

        # Verify that the first gene document looks as expected
        mock_dir = 'matrix_mtx'
        doc, expected_doc = get_nth_gene_docs(0, docs, mock_dir)
        self.assertEqual(doc, expected_doc)

    def test_mtx_bundle_argument_validation(self):
        """Omitting --gene-file and --barcode-file in MTX ingest should error
        """

        args = [
            '--study-id',
            'SCP1',
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/matrix.mtx',
            '--matrix-file-type',
            'mtx',
        ]

        self.assertRaises(ValueError, self.setup_ingest, args)

    def test_ingest_loom(self):
        """Ingest Pipeline should extract, transform, and load Loom files
        """

        study_accession = get_random_study_accession()

        args = [
            '--study-id',
            study_accession,
            '--study-file-id',
            '1234abc',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/test_loom.loom',
            '--matrix-file-type',
            'loom',
        ]

        ingest = self.setup_ingest(args)

        genes = ingest.db.collection(u'genes')
        stream = (
            genes.where(u'taxon_common_name', u'==', u'human')
            .where(u'study_accession', u'==', study_accession)
            .stream()
        )

        docs = [doc.to_dict() for doc in stream]

        # Verify that 10 gene docs were written to Firestore
        num_docs = len(docs)
        self.assertEqual(num_docs, 10)

        # Verify that the first gene document looks as expected
        mock_dir = 'loom'
        doc, expected_doc = get_nth_gene_docs(0, docs, mock_dir)
        self.assertEqual(doc, expected_doc)


if __name__ == '__main__':
    unittest.main()
