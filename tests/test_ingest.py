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

class IngestTestCase(unittest.TestCase):

    def test_ingest_dense(self):
        """Ingest Service should handle dense matrices
        """
        output_dir = 'test_output/'

        args = [
            'ingest_expression',
            '--matrix-file', '../tests/data/dense_matrix_19_genes_100k_cells.txt',
            '--matrix-file-type', 'dense'
        ]
        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        IngestService.load_expression_data = mock_load_expression_data

        mock_db = MockFirestore()

        ingest = IngestService(**arguments, db=mock_db)

        if hasattr(ingest, 'ingest_expression'):
            getattr(ingest, 'ingest_expression')()
        
        expression_models = ingest.load_expression_data_args[0]

        num_models = len(expression_models)
        expected_num_models = 19
        self.assertEqual(num_models, expected_num_models)

        first_model = expression_models[0].__dict__
        with open('mock_data/dense_matrix_19_genes_100k_cells/gene_model_1.json') as f:
            expected_first_model = json.loads(f.read())

        self.assertEqual(first_model, expected_first_model)

if __name__ == '__main__':
    unittest.main()