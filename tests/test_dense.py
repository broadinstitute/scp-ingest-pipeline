import unittest
import sys

sys.path.append("../ingest")
from dense import Dense


class TestDense(unittest.TestCase):
    def test_validate_gene_keyword_false(self):
        "Validates validate_gene_keyword() returns false correctly"
        expression_matrix = Dense(
            '../tests/data/expression_matrix_bad_missing_keyword.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertFalse(expression_matrix.validate_gene_keyword())
        expression_matrix.close()

    def test_validate_gene_keyword_true(self):
        "Validates validate_gene_keyword() returns true correctly"
        expression_matrix = Dense(
            '../tests/data/dense_expression_matrix.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertTrue(expression_matrix.validate_gene_keyword())
        expression_matrix.close()

    def test_validate_unique_header_false(self):
        "Validates validate_unique_header() returns false correctly"
        expression_matrix = Dense(
            '../tests/data/expression_matrix_non_unique.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertFalse(expression_matrix.validate_unique_header())
        expression_matrix.close()

    def test_validate_unique_header_true(self):
        "Validates validate_unique_header() returns true correctly"
        expression_matrix = Dense(
            '../tests/data/dense_expression_matrix.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertTrue(expression_matrix.validate_unique_header())
        expression_matrix.close()

    def test_duplicate_gene(self):
        expression_matrix = Dense(
            '../tests/data/expression_matrix_duplicate_header.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        expression_matrix.preprocess()
        with self.assertRaises(ValueError):
            for gene in expression_matrix.transform():
                pass
        expression_matrix.close()

    def test_validate_format_false(self):
        expression_matrix = Dense(
            '../tests/data/expression_matrix_bad_missing_keyword.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertFalse(expression_matrix.validate_format())
        expression_matrix.close()

    def test_validate_format_true(self):
        expression_matrix = Dense(
            '../tests/data/dense_expression_matrix.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            tracer=None,
        )
        self.assertTrue(expression_matrix.validate_format())
