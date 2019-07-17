"""Unit tests for ingest.py

EXAMPLES
cd tests; python3 test_ingest.py
"""

import unittest
from glob import glob
import json
import sys

from mockfirestore import MockFirestore

sys.path.append('../ingest')
from ingest_pipeline import *

def mock_load_expression_data(self, *args, **kwargs):
    self.load_expression_data_args = args
    self.load_expression_data_kwargs = kwargs


IngestService.load_expression_data = mock_load_expression_data


class IngestTestCase(unittest.TestCase):

    def setup_ingest(self, args):
        args_list = args.split(' ')
        parsed_args = create_parser().parse_args(args_list)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        mock_db = MockFirestore()

        ingest = IngestService(**arguments, db=mock_db)

        if hasattr(ingest, 'ingest_expression'):
            getattr(ingest, 'ingest_expression')()

        return ingest

    def test_ingest_dense_matrix(self):
        """Ingest Service should handle dense matrices
        """

        args = ('ingest_expression '
                '--matrix-file ../tests/data/dense_matrix_19_genes_100k_cells.txt '
                '--matrix-file-type dense')
        ingest = self.setup_ingest(args)

        expression_models = ingest.load_expression_data_args[0]

        num_models = len(expression_models)
        expected_num_models = 19
        self.assertEqual(num_models, expected_num_models)

        first_model = expression_models[0].__dict__
        with open('mock_data/dense_matrix_19_genes_100k_cells_txt/gene_model_1.json') as f:
            expected_first_model = json.loads(f.read())

        self.assertEqual(first_model, expected_first_model)

    def test_ingest_mtx_matrix(self):
        """Ingest Service should handle MTX matrix bundles
        """

        args = ('ingest_expression '
                '--matrix-file ../tests/data/matrix.mtx '
                '--matrix-file-type mtx '
                '--gene-file ../tests/data/genes.tsv '
                '--barcode-file ../tests/data/barcodes.tsv')
        ingest = self.setup_ingest(args)

        expression_models = ingest.load_expression_data_args[0]

        num_models = len(expression_models)
        expected_num_models = 25
        self.assertEqual(num_models, expected_num_models)

        first_model = list(expression_models)[0].__dict__
        with open('mock_data/matrix_mtx/gene_model_1.json') as f:
            expected_first_model = json.loads(f.read())
        self.assertEqual(first_model, expected_first_model)


if __name__ == '__main__':
    unittest.main()
