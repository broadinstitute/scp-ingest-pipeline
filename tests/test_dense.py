import sys
import unittest
from unittest.mock import patch, MagicMock
from bson.objectid import ObjectId

sys.path.append("../ingest")
from expression_files.dense_ingestor import DenseIngestor
from mock_data.dense_matrix_19_genes_100k_cells_txt.gene_models_0 import gene_models


class TestDense(unittest.TestCase):
    def test_process_header(self):
        header = ["one", "two", '"three', "'four"]
        actual_header = DenseIngestor.process_header(header)
        self.assertEqual(["one", "two", "three", "four"], actual_header)

    def test_check_header_valid_values(self):
        header = ["one", "two", '"three', "'four"]
        self.assertTrue(DenseIngestor.check_header_valid_values(header))

        empty_value = ["", "two"]
        self.assertRaises(
            ValueError, DenseIngestor.check_header_valid_values, empty_value
        )
        space_value = ["   ", "two"]
        self.assertRaises(
            ValueError, DenseIngestor.check_header_valid_values, space_value
        )
        nan_value = ["nan", "two"]
        self.assertRaises(
            ValueError, DenseIngestor.check_header_valid_values, nan_value
        )

    def test_filter_expression_scores(self):
        scores = ["BRCA1", 4, 0, 3, "0", None, "", "   ", "nan", "NaN", "Nan", "NAN"]
        cells = ["foo", "foo2", "foo3", "foo4", "foo5", "foo6", "foo7"]
        actual_filtered_values, actual_filtered_cells = DenseIngestor.filter_expression_scores(
            scores[1:], cells
        )
        self.assertEqual([4, 3], actual_filtered_values)
        self.assertEqual(["foo2", "foo4"], actual_filtered_cells)

        invalid_scores = ["BRCA1", 4, 0, 3, "0", None, "T"]
        self.assertRaises(
            ValueError, DenseIngestor.filter_expression_scores, invalid_scores, cells
        )

    def test_process_row(self):
        # Positive test case
        valid_row = ["BRCA1", "' 1.45678 '", '"3.45678"', "2"]
        processed_row = DenseIngestor.process_row(valid_row)
        self.assertEqual([1.457, 3.457, 2.0], processed_row)

        # Negative test case
        invalid_row = ["BRCA1", "' 1.BRCA1 '", '"3.45678"', "2"]
        self.assertRaises(ValueError, DenseIngestor.process_row, invalid_row)

    def test_check_gene_keyword(self):
        """Validates validate_gene_keyword() returns false correctly"""
        # Mimics the row following the header
        row = ["BRCA1", "' 1.45678 '", '"3.45678"', "2"]

        header = ["GENE", "foo", "foo2", "foo3"]
        r_header = ["foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.check_gene_keyword(header, row))
        self.assertTrue(DenseIngestor.check_gene_keyword(r_header, row))

        empty_header = [""]
        invalid_header = ["foo", "foo2", "foo3", "foo4"]
        self.assertRaises(ValueError, DenseIngestor.check_gene_keyword, invalid_header, row)
        self.assertRaises(ValueError, DenseIngestor.check_gene_keyword, empty_header, row)

    def test_check_unique_header(self):
        """Validates validate_unique_header() returns false correctly"""
        header = ["GENE", "foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.check_unique_header(header))

        invalid_header = ["GENE", "foo", "foo", "foo3"]
        self.assertRaises(ValueError, DenseIngestor.check_unique_header, invalid_header)

    def test_duplicate_gene(self):
        expression_matrix = DenseIngestor(
            "../tests/data/expression_matrix_bad_duplicate_gene.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            tracer=None,
        )
        with self.assertRaises(ValueError) as error:
            for gene_mode, data_array in expression_matrix.transform():
                pass
        self.assertEqual(str(error.exception), "Duplicate gene: Itm2a")

    @patch("expression_files.expression_files.GeneExpression.check_unique_cells")
    def test_check_valid(
        self,
        mock_check_unique_cells
    ):
        """Confirms the errors are correctly promulgated"""
        query_params = (ObjectId("5dd5ae25421aa910a723a337"), MagicMock())

        self.assertTrue(
            DenseIngestor.check_valid(
                ["GENE", "foo", "foo2"],
                ["gene1", "0", "1"],
                query_params
            )
        )

        with self.assertRaises(ValueError) as cm:
            DenseIngestor.check_valid(
                ["GENE", "foo", "foo"],
                ["gene1", "0", "1"],
                query_params
            )
        expected_msg = "Duplicate header values are not allowed"
        self.assertEqual(expected_msg, str(cm.exception))

        # Should concatenate the two value errors
        mock_check_unique_cells.return_value = True
        with self.assertRaises(ValueError) as cm:
            DenseIngestor.check_valid(
                ["foo", "nan"],
                ["foo2", "foo3"],
                query_params
            )
        expected_msg = "Required 'GENE' header is not present; nan is not allowed as a header value"
        self.assertEqual(expected_msg, str(cm.exception))

    @patch("expression_files.expression_files.GeneExpression.load")
    @patch(
        "expression_files.dense_ingestor.DenseIngestor.transform",
        return_value=[("foo1", "foo2")],
    )
    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_execute_ingest(self, mock_load, mock_transform, mock_has_unique_cells):
        """
        Integration test for execute_ingest()
        """
        expression_matrix = DenseIngestor(
            "../tests/data/expression_matrix_bad_duplicate_gene.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        # When is_valid_format() is false exception should be raised
        self.assertRaises(ValueError, expression_matrix.execute_ingest())

        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        expression_matrix.execute_ingest()
        self.assertTrue(mock_transform.called)
        self.assertTrue(mock_load.called)

    def test_transform_fn(self):
        """
        Assures transform function creates gene data model correctly
        """
        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        for actual_gene_models, actual_data_arrays in expression_matrix.transform():
            # _id is a unique identifier and can not be predicted
            # so we exclude it from the comparison
            for actual_gene_model in actual_gene_models:
                del actual_gene_model["_id"]
                gene_name = actual_gene_model["name"]
                self.assertEqual(actual_gene_model, gene_models[gene_name])
