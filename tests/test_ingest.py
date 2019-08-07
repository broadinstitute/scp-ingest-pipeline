"""Unit tests for ingest.py

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
from glob import glob
import sys
import os
import unittest
from unittest.mock import patch
from shutil import copyfile

import google.cloud

sys.path.append('../ingest')
from ingest_pipeline import *

def mock_load_expression_data(self, *args, **kwargs):
    """Enables overwriting normal function with this placeholder.
    Returning the arguments enables tests to verify that the code invokes
    this method with expected argument values.

    Ideally we would mock Firestore itself, but existing libraries to mock
    Firestore (e.g. https://github.com/mdowds/python-mock-firestore) don't
    implement important Firestore methods like "batch()".
    """
    self.load_expression_data_args = args
    self.load_expression_data_kwargs = kwargs

def mock_storage_client():

    class MockStorageBucket():
        def __init__(self, name):
            self.name = name
            return

        def blob(self, blob_name):
            return mock_storage_blob(bucket=self.name, name=blob_name)

    class MockStorageClient():
        def __init__(self):
            return

        def get_bucket(bucket_name):
            return MockStorageBucket(bucket_name)

    return MockStorageClient

def mock_storage_blob(*args, **kwargs):
    """Mocks Google Cloud Storage library

    TODO: Watch progress on official Storage emulator for integration tests:
        - https://github.com/googleapis/google-cloud-python/issues/8728
        - https://github.com/googleapis/google-cloud-python/issues/4840

    When such an emulator is released, use it and remove this custom mock.
    """

    class MockStorageBlob():
        def __init__(self, bucket=None, name=None):
            self.bucket = bucket
            self.name = '../' + name

        def exists(self, storage_client):
            return os.path.exists(self.name)

        def download_to_filename(self, filename):
            """Mock; doesn't actually download.  Makes local copy instead."""
            copyfile(self.name, filename)

    return MockStorageBlob(*args, **kwargs)

def mock_firestore_client():
    """Mocks firestore.Client() by returning nothing upon initializing client

    See notes in mock_load_expression_data for context.
    """
    return

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
    # print('actual_model')
    # print(actual_model)

    with open(f'mock_data/{mock_dir}/gene_model_{n}.txt') as f:
        # Create a dictionary from the string-literal mock
        expected_model = ast.literal_eval(f.read())

    return actual_model, expected_model

# Mock method that loads data to Firestore
IngestPipeline.load_expression_data = mock_load_expression_data

class IngestTestCase(unittest.TestCase):

    @patch('google.cloud.storage.Blob', side_effect=mock_storage_blob)
    @patch('google.cloud.storage.Client', side_effect=mock_storage_client)
    @patch('google.cloud.firestore.Client', side_effect=mock_firestore_client)
    def setup_ingest(self, args, mock_firestore_client, mock_storage_client, mock_storage_blob):
        args_list = args.split(' ')
        parsed_args = create_parser().parse_args(args_list)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        ingest = IngestPipeline(**arguments)

        if hasattr(ingest, 'ingest_expression'):
            getattr(ingest, 'ingest_expression')()

        return ingest

    def test_ingest_dense_matrix(self):
        """Ingest Pipeline should extract and transform remote file
        """

        args = ('ingest_expression '
                '--matrix-file gs://fake-bucket/tests/data/dense_matrix_19_genes_100k_cells.txt '
                '--matrix-file-type dense')
        ingest = self.setup_ingest(args)

        models = ingest.load_expression_data_args[0]

        # Verify that 19 gene models were passed into load method
        num_models = len(models)
        expected_num_models = 19
        self.assertEqual(num_models, expected_num_models)

        # Verify that the first gene model looks as expected
        mock_dir = 'dense_matrix_19_genes_100k_cells_txt'
        model, expected_model = get_nth_gene_models(0, models, mock_dir)

        # Adjust for change needed per mock
        expected_source_file_name =\
            'gs://fake-bucket/tests/data/dense_matrix_19_genes_100k_cells.txt'
        actual_source_file_name = model['source_file_name'][0]  # Why does app code make this a tuple?
        self.assertEqual(actual_source_file_name, expected_source_file_name)

        model['source_file_name'] = ('/tests/data/dense_matrix_19_genes_100k_cells.txt',)  # Why does app code make this a tuple?
        model['subdocument']['source_file_name'] = '/tests/data/dense_matrix_19_genes_100k_cells.txt'

        self.assertEqual(model, expected_model)

    def test_ingest_local_dense_matrix(self):
        """Ingest Pipeline should extract and transform dense matrices
        """

        args = ('ingest_expression '
                '--matrix-file ../tests/data/dense_matrix_19_genes_100k_cells.txt '
                '--matrix-file-type dense')
        ingest = self.setup_ingest(args)

        models = ingest.load_expression_data_args[0]

        # Verify that 19 gene models were passed into load method
        num_models = len(models)
        expected_num_models = 19
        self.assertEqual(num_models, expected_num_models)

        # Verify that the first gene model looks as expected
        mock_dir = 'dense_matrix_19_genes_100k_cells_txt'
        model, expected_model = get_nth_gene_models(0, models, mock_dir)

        self.assertEqual(model, expected_model)

    def test_ingest_mtx_matrix(self):
        """Ingest Pipeline should extract and transform MTX matrix bundles
        """

        args = ('ingest_expression '
                '--matrix-file ../tests/data/matrix.mtx '
                '--matrix-file-type mtx '
                '--gene-file ../tests/data/genes.tsv '
                '--barcode-file ../tests/data/barcodes.tsv')
        ingest = self.setup_ingest(args)

        models = ingest.load_expression_data_args[0]

        # Verify that 25 gene models were passed into load method
        num_models = len(models)
        expected_num_models = 25
        self.assertEqual(num_models, expected_num_models)

        # Verify that the first gene model looks as expected
        mock_dir = 'matrix_mtx'
        model, expected_model = get_nth_gene_models(0, models, mock_dir)
        self.assertEqual(model, expected_model)

    def test_mtx_bundle_argument_validation(self):
        """Omitting --gene-file and --barcode-file in MTX ingest should error
        """

        args = ('ingest_expression '
                '--matrix-file ../tests/data/matrix.mtx '
                '--matrix-file-type mtx')


        self.assertRaises(ValueError, self.setup_ingest, args)

if __name__ == '__main__':
    unittest.main()
