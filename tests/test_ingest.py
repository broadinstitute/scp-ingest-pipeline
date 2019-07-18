"""Unit tests for ingest.py

EXAMPLES
cd tests; python3 test_ingest.py
"""

from glob import glob
import json
import sys
import unittest
from unittest.mock import patch

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

def mock_firestore_client(*args, **kwargs):
    """Return nothing upon invoking firestore.Client()
    """
    return

def get_nth_gene_models(n, models, mock_dir):
    """Return Nth actual and expected gene models, using actual and mock data
    """
    actual_model = expression_models[0].__dict__
    path = f'mock_data/${mock_path}/gene_model_${n}.json'
    with open(path) as f:
        expected_model = json.loads(f.read())
    return actual_model, expected_model

IngestPipeline.load_expression_data = mock_load_expression_data

class IngestTestCase(unittest.TestCase):

    def setup_ingest(self, args):
        args_list = args.split(' ')
        parsed_args = create_parser().parse_args(args_list)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        ingest = IngestPipeline(**arguments)

        if hasattr(ingest, 'ingest_expression'):
            getattr(ingest, 'ingest_expression')()

        return ingest

    @patch('google.cloud.firestore.Client', side_effect=mock_firestore_client)
    def test_ingest_dense_matrix(self, mock_firestore_client):
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

    @patch('google.cloud.firestore.Client', side_effect=mock_firestore_client)
    def test_ingest_mtx_matrix(self, mock_firestore_client):
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


if __name__ == '__main__':
    unittest.main()
