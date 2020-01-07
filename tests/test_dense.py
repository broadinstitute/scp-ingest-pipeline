import unittest
import sys

sys.path.append("../ingest")
from dense import Dense


class TestAnnotations(unittest.TestCase):
    def test_validate_gene_keyword(self):
        expression_matrix = Dense(
            '../tests/data/expression_matrix_bad_missing_keyword.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        with self.assertRaises(ValueError):
            expression_matrix.validate_gene_keyword()

    def test_unique_header(self):
        expression_matrix = Dense(
            '../tests/data/expression_matrix_non-unique.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        with self.assertRaises(ValueError):
            expression_matrix.validate_unique_header()

    def test_duplicate_gene(self):
        expression_matrix = Dense(
            '../tests/data/expression_matrix_duplicate_header.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        with self.assertRaises(ValueError):
            for gene in expression_matrix.transform():
                pass
