import sys
import unittest
from unittest.mock import patch

sys.path.append("../ingest")
from expression_files.dense_ingestor import DenseIngestor
from mock_data.dense_matrix_19_genes_100k_cells_txt.gene_models_0 import gene_models
from mock_data.dense_matrix_19_genes_100k_cells_txt.data_arrays import data_arrays


class TestDense(unittest.TestCase):
    # TODO add test method after SCP-2635 implemented

    def test_process_header(self):
        header = ["one", "two", '"three', "'four"]
        actual_header = DenseIngestor.process_header(header)
        self.assertEqual(["one", "two", "three", "four"], actual_header)

    def test_header_has_valid_values(self):
        header = ["one", "two", '"three', "'four"]
        self.assertTrue(DenseIngestor.header_has_valid_values(header))

        empty_value = ["", "two"]
        space_value = ["   ", "two"]
        nan_value = ["nan", "two"]
        self.assertFalse(DenseIngestor.header_has_valid_values(empty_value))
        self.assertFalse(DenseIngestor.header_has_valid_values(space_value))
        self.assertFalse(DenseIngestor.header_has_valid_values(nan_value))

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
        valid_row = ["' 1.45678 '", '"3.45678"', "2"]
        processed_row = DenseIngestor.process_row(valid_row)
        self.assertEqual([1.457, 3.457, 2.0], processed_row)

        # Negative test case
        invalid_row = ["' 1.BRCA1 '", '"3.45678"', "2"]
        self.assertRaises(ValueError, DenseIngestor.process_row, invalid_row)

    def test_has_gene_keyword(self):
        """Validates validate_gene_keyword() returns false correctly"""
        # Mimics the row following the header
        row = ["BRCA1", "' 1.45678 '", '"3.45678"', "2"]

        header = ["GENE", "foo", "foo2", "foo3"]
        r_header = ["foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.has_gene_keyword(header, row))
        self.assertTrue(DenseIngestor.has_gene_keyword(r_header, row))

        empty_header = [""]
        invalid_header = ["foo", "foo2", "foo3", "foo4"]
        self.assertFalse(DenseIngestor.has_gene_keyword(invalid_header, row))
        self.assertFalse(DenseIngestor.has_gene_keyword(empty_header, row))

    def test_has_unique_header(self):
        """Validates validate_unique_header() returns false correctly"""
        header = ["GENE", "foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.has_unique_header(header))

        invalid_header = ["GENE", "foo", "foo", "foo3"]
        self.assertFalse(DenseIngestor.has_unique_header(invalid_header))

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

    @patch("expression_files.dense_ingestor.DenseIngestor.has_unique_header")
    @patch("expression_files.dense_ingestor.DenseIngestor.has_gene_keyword")
    @patch("expression_files.dense_ingestor.DenseIngestor.header_has_valid_values")
    def test_is_valid_format(
        self,
        mock_has_unique_header,
        mock_has_gene_keyword,
        mock_header_has_valid_values,
    ):
        """Confirms functions in is_valid_format() are called"""
        # Should raise Value error
        mock_has_unique_header.return_value = False
        mock_has_gene_keyword.return_value = False
        mock_header_has_valid_values.return_value = False
        self.assertFalse(
            DenseIngestor.is_valid_format(["foo", "foo1"], ["foo2", "foo3"])
        )
        self.assertTrue(mock_has_unique_header.called)
        self.assertTrue(mock_has_gene_keyword.called)
        self.assertTrue(mock_header_has_valid_values.called)

        # Should raise Value error
        mock_has_unique_header.return_value = True
        mock_has_gene_keyword.return_value = False
        mock_header_has_valid_values.return_value = False
        self.assertFalse(
            DenseIngestor.is_valid_format(["foo", "foo1"], ["foo2", "foo3"])
        )

        # Should raise Value error
        mock_has_unique_header.return_value = True
        mock_has_gene_keyword.return_value = True
        mock_header_has_valid_values.return_value = False
        self.assertFalse(
            DenseIngestor.is_valid_format(["foo", "foo1"], ["foo2", "foo3"])
        )

        # When is_valid_format() returns true
        mock_has_unique_header.return_value = True
        mock_has_gene_keyword.return_value = True
        mock_header_has_valid_values.return_value = True
        self.assertTrue(
            DenseIngestor.is_valid_format(["foo", "foo1"], ["foo2", "foo3"])
        )

    @patch("expression_files.expression_files.GeneExpression.load")
    @patch(
        "expression_files.dense_ingestor.DenseIngestor.transform",
        return_value=[("foo1", "foo2")],
    )
    def test_execute_ingest(self, mock_load, mock_transform):
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
            for actual_data_array in actual_data_arrays:
                del actual_data_array["linear_data_id"]
                data_array_name = actual_data_array["name"]
                self.assertEqual(actual_data_array, data_arrays[data_array_name])
