import sys
import unittest
from unittest.mock import patch, MagicMock, PropertyMock
from bson.objectid import ObjectId

sys.path.append("../ingest")
from expression_files.dense_ingestor import DenseIngestor
from expression_files.expression_files import GeneExpression
from mock_data.dense_matrix_19_genes_100k_cells_txt.gene_models_0 import gene_models
from mock_data.dense_matrix_19_genes_100k_cells_txt.data_arrays import data_arrays


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
        scores = ["BRCA1", 4, 0, "0", None, "", "   ", "nan", "NaN", "Nan", "NAN", 3]
        cells = [
            "foo",
            "foo2",
            "foo3",
            "foo4",
            "foo5",
            "foo6",
            "foo7",
            "foo8",
            "foo9",
            "foo10",
            "foo11",
        ]
        actual_filtered_values, actual_filtered_cells = DenseIngestor.filter_expression_scores(
            scores[1:], cells
        )
        self.assertEqual([4, 3], actual_filtered_values)
        self.assertEqual(["foo", "foo11"], actual_filtered_cells)

        invalid_scores = ["BRCA1", 4, 0, 3, "0", None, "T"]
        self.assertRaises(
            ValueError, DenseIngestor.filter_expression_scores, invalid_scores, cells
        )
        not_enough_scores = ["BRCA1", 4, 0, "0"]
        with self.assertRaises(ValueError) as error:
            DenseIngestor.filter_expression_scores(not_enough_scores, cells)

        too_many_scores = ["BRCA1", 4, 0, "0", 5, 6]
        with self.assertRaises(ValueError) as error:
            DenseIngestor.filter_expression_scores(too_many_scores, cells)
        expected_msg = (
            "Number of cell and expression values must be the same. "
            f"Found row with {len(too_many_scores)} expression values."
            f"Header contains {len(cells)} cells"
        )
        self.assertEqual(expected_msg, str(error.exception))

    def test_process_row(self):
        # Positive test case
        valid_row = ["' 1.45678 '", '"3.45678"', "2"]
        processed_row = DenseIngestor.process_row(valid_row)
        self.assertEqual([1.457, 3.457, 2.0], processed_row)

        # Negative test case
        invalid_row = ["' 1.BRCA1 '", '"3.45678"', "2"]
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
        self.assertRaises(
            ValueError, DenseIngestor.check_gene_keyword, invalid_header, row
        )
        self.assertRaises(
            ValueError, DenseIngestor.check_gene_keyword, empty_header, row
        )

    def test_is_r_formatted_file(self):
        # Mimics the row following the header
        row = ["BRCA1", "' 1.45678 '", '"3.45678"', "2"]

        r_header = ["foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.is_r_formatted_file(r_header, row))

        dense_header = ["GENE", "foo", "foo2", "foo3"]
        self.assertFalse(DenseIngestor.is_r_formatted_file(dense_header, row))

    def test_check_unique_header(self):
        """Validates validate_unique_header() returns false correctly"""
        header = ["GENE", "foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.check_unique_header(header))

        invalid_header = ["GENE", "foo", "foo", "foo3"]
        self.assertRaises(ValueError, DenseIngestor.check_unique_header, invalid_header)

    def test_duplicate_gene(self):
        FILE_PATH = "../tests/data/expression_matrix_bad_duplicate_gene.txt"
        expression_matrix = DenseIngestor(
            FILE_PATH,
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            tracer=None,
        )
        # Return file handler to correct position
        expression_matrix.csv_file_handler = expression_matrix.open_file(FILE_PATH)[0]
        # Skip header
        next(expression_matrix.csv_file_handler)
        with self.assertRaises(ValueError) as error:
            for gene_mode, data_array in expression_matrix.transform():
                pass
        self.assertEqual(str(error.exception), "Duplicate gene: Itm2a")

    @patch("expression_files.expression_files.GeneExpression.check_unique_cells")
    def test_check_valid(self, mock_check_unique_cells):
        """Confirms the errors are correctly promulgated"""
        query_params = (ObjectId("5dd5ae25421aa910a723a337"), MagicMock())

        self.assertTrue(
            DenseIngestor.check_valid(
                ["GENE", "foo", "foo2"], ["gene1", "0", "1"], query_params
            )
        )

        with self.assertRaises(ValueError) as cm:
            DenseIngestor.check_valid(
                ["GENE", "foo", "foo"], ["gene1", "0", "1"], query_params
            )
        expected_msg = "Duplicate header values are not allowed"
        self.assertEqual(expected_msg, str(cm.exception))

        # Should concatenate the two value errors
        mock_check_unique_cells.return_value = True
        with self.assertRaises(ValueError) as cm:
            DenseIngestor.check_valid(["foo", "nan"], ["foo2", "foo3"], query_params)
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
            for actual_data_array in actual_data_arrays:
                del actual_data_array["linear_data_id"]
                data_array_name = actual_data_array["name"]
                self.assertEqual(actual_data_array, data_arrays[data_array_name])

    def test_transform_fn_batch(self):
        """
        Assures transform function batches data array creation
        the +1 fudge factor is because we only check for batch size after a full row
        has been processed
        """

        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_10_genes_15_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        with patch(
            "expression_files.expression_files.GeneExpression.DATA_ARRAY_BATCH_SIZE",
            new_callable=PropertyMock,
            return_value=4,
        ):
            for actual_gene_models, actual_data_arrays in expression_matrix.transform():
                self.assertTrue(
                    GeneExpression.DATA_ARRAY_BATCH_SIZE + 1 >= len(actual_data_arrays)
                )

        """
        Assures transform function creates gene data model correctly for files
        with a number of data arrays that is a multiple of DATA_ARRAY_BATCH_SIZE (SCP-2669)
        (note that the number of data_arrays is (number of genes) * 2 + 1)
        """
        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_10_genes_15_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        with patch(
            "expression_files.expression_files.GeneExpression.DATA_ARRAY_BATCH_SIZE",
            new_callable=PropertyMock,
            return_value=21,
        ):
            for actual_gene_models, actual_data_arrays in expression_matrix.transform():
                self.assertEqual(
                    GeneExpression.DATA_ARRAY_BATCH_SIZE, len(actual_data_arrays)
                )
