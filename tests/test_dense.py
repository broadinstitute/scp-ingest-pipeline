import sys
import unittest
from unittest.mock import patch, MagicMock, PropertyMock
from bson.objectid import ObjectId
import csv

from test_expression_files import mock_load_genes_batched, mock_expression_load

# R models
from mock_data.expression.r_format.models import r_data_models

sys.path.append("../ingest")
from expression_files.dense_ingestor import DenseIngestor
from expression_files.expression_files import GeneExpression
from ingest_files import DataArray


def mock_load_no_exp_data(documents, collection_name):
    """
    Overwrites GeneExpression.load()
    Confirms gene models are created when a gene does not have expression values
    """
    if collection_name == DataArray.COLLECTION_NAME:
        # There will always be a data array model for the cell names
        assert len(documents) == 1
    else:
        assert collection_name == GeneExpression.COLLECTION_NAME
        assert len(documents) == 4


def mock_load_r_files(documents, collection_name):
    """
    Overwrites GeneExpression.load() for R formatted file

    GeneExpression.load() is called multiple times. This method will verify
    models in the arguments have the expected values.
    """
    # _id and linear_data_id are unique identifiers and can not be predicted
    # so we exclude it from the comparison
    for document in documents:
        model_name = document["name"]
        if collection_name == GeneExpression.COLLECTION_NAME:
            # Ensure that 'ObjectID' in model is removed
            del document["_id"]
            # Verify gene model looks as expected
            assert document == r_data_models["gene_models"][model_name]
        else:
            del document["linear_data_id"]
            assert document == r_data_models["data_arrays"][model_name]


class TestDense(unittest.TestCase):
    GeneExpression.load = mock_expression_load

    def test_set_header(self):
        # R File where last value is ""
        with open("../tests/data/r_format_text.txt") as f:
            r_file = csv.reader(f, delimiter="\t")
            header = DenseIngestor.set_header(r_file)
            # Reset file handle to grab expected header
            f.seek(0)
            expected_header = next(r_file)[:-1]
            self.assertEqual(header, expected_header)

        # R File where "GENE" isn't present
        with open("../tests/data/r_format_text.csv") as f:
            r_file = csv.reader(f)
            header = DenseIngestor.set_header(r_file)
            # Reset file handle to grab expected header
            f.seek(0)
            expected_header = next(r_file)
            self.assertEqual(header, expected_header)

        # Regular dense matrix
        with open("../tests/data/dense_matrix_10_genes_15_cells.txt") as f:
            dense_file = csv.reader(f)
            header = DenseIngestor.set_header(dense_file)
            # Reset file handle to grab expected header
            f.seek(0)
            expected_header = next(dense_file)[1:]
            self.assertEqual(header, expected_header)

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
        (
            actual_filtered_values,
            actual_filtered_cells,
        ) = DenseIngestor.filter_expression_scores(scores[1:], cells)
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
            f"Found row with {len(too_many_scores)} expression values. "
            f"Header contains {len(cells)} cells."
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
        """Validates validate_gene_keyword() returns false correctly."""
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
        is_r, header = DenseIngestor.is_r_formatted_file(r_header, row)
        self.assertTrue(is_r)
        self.assertEqual(header, r_header)

        r_header = ["foo", "foo2", "foo3", ""]
        is_r, header = DenseIngestor.is_r_formatted_file(r_header, row)
        self.assertTrue(is_r)
        self.assertEqual(header, ["foo", "foo2", "foo3"])

        dense_header = ["GENE", "foo", "foo2", "foo3"]
        is_r, header = DenseIngestor.is_r_formatted_file(dense_header, row)
        self.assertFalse(is_r)
        self.assertEqual(header, ["foo", "foo2", "foo3"])

    def test_check_unique_header(self):
        """Validates validate_unique_header() returns false correctly."""
        header = ["GENE", "foo", "foo2", "foo3"]
        self.assertTrue(DenseIngestor.check_unique_header(header))

        invalid_header = ["GENE", "foo", "foo", "foo3"]
        self.assertRaises(ValueError, DenseIngestor.check_unique_header, invalid_header)

    @patch(
        "expression_files.expression_files.GeneExpression.is_raw_count_file",
        return_value=False,
    )
    def test_duplicate_gene(self, moock_is_raw_count_file):
        FILE_PATH = "../tests/data/expression_matrix_bad_duplicate_gene.txt"
        expression_matrix = DenseIngestor(
            FILE_PATH,
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            tracer=None,
        )
        with self.assertRaises(ValueError) as error:
            for gene_mode, data_array in expression_matrix.transform():
                pass
        self.assertEqual(str(error.exception), "Duplicate gene: Itm2a")

    @patch("expression_files.expression_files.GeneExpression.check_unique_cells")
    def test_check_valid(self, mock_check_unique_cells):
        """Confirms the errors are correctly promulgated."""
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
        expected_msg = (
            'Duplicate header values are not allowed. Duplicates include: "foo"'
        )
        self.assertEqual(expected_msg, str(cm.exception))

        # Should concatenate the two value errors
        mock_check_unique_cells.return_value = True
        with self.assertRaises(ValueError) as cm:
            DenseIngestor.check_valid(["foo", "nan"], ["foo2", "foo3"], query_params)
        expected_msg = 'Required "GENE" header is not present.; "nan" is not allowed as a header value.'
        self.assertEqual(expected_msg, str(cm.exception))

    # Commenting out "bad_execute_ingest" validation due to change in mock with
    # Python >= 3.8 (accessing attributes on the mock object before patching
    # mock object, causing error not seen with Python 3.7)
    # TODO (SCP-5032): Resolve this test
    #
    # @patch(
    #     "expression_files.expression_files.GeneExpression.is_raw_count_file",
    #     return_value=False,
    # )
    # @patch(
    #     "expression_files.expression_files.GeneExpression.check_unique_cells",
    #     return_value=True,
    # )
    # @patch(
    #     "expression_files.dense_ingestor.DenseIngestor.transform",
    #     return_value=[("foo1", "foo2")],
    # )
    # @patch("expression_files.expression_files.GeneExpression.load")
    # def test_bad_execute_ingest(
    #     self, mock_load, mock_transform, mock_has_unique_cells, mock_is_raw_count_file
    # ):
    #     """
    #     Integration test for bad execute_ingest()
    #     """
    #     expression_matrix = DenseIngestor(
    #         "../tests/data/expression_matrix_bad_duplicate_gene.txt",
    #         "5d276a50421aa9117c982844",
    #         "5dd5ae25421aa910a723a336",
    #     )
    #     # When is_valid_format() is false exception should be raised
    #     self.assertRaises(ValueError, expression_matrix.execute_ingest())

    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    @patch(
        "expression_files.dense_ingestor.DenseIngestor.transform",
        return_value=[("foo1", "foo2")],
    )
    @patch("expression_files.expression_files.GeneExpression.load")
    def test_good_execute_ingest(
        self, mock_load, mock_transform, mock_has_unique_cells
    ):
        """
        Integration test for good execute_ingest()
        """
        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        expression_matrix.execute_ingest()
        self.assertTrue(mock_transform.called)
        # Commenting out this "good_execute_ingest" assertion due to change in
        # mock with Python >= 3.8 (accessing attributes on the mock object
        # before patching mock object, causing error not seen with Python 3.7)
        # TODO (SCP-5032): Resolve this test
        #
        # self.assertTrue(mock_load.called)

    @patch("expression_files.expression_files.GeneExpression.is_raw_count_file")
    def test_transform_fn(self, mock_is_raw_count_file):
        """
        Assures transform function creates data models correctly.
        """
        mock_is_raw_count_file.return_value = False
        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        expression_matrix.test_models = None
        expression_matrix.models_processed = 0
        expression_matrix.transform()
        amount_of_models = len(
            expression_matrix.test_models["data_arrays"].keys()
        ) + len(expression_matrix.test_models["gene_models"].keys())
        self.assertEqual(expression_matrix.models_processed, amount_of_models)

        # Test raw count files
        mock_is_raw_count_file.return_value = True
        expression_matrix = DenseIngestor(
            "../tests/data/raw1_human_5k_cells_80_genes.dense.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        expression_matrix.test_models = None
        expression_matrix.models_processed = 0
        expression_matrix.transform()
        amount_of_models = len(expression_matrix.test_models["data_arrays"].keys())
        self.assertEqual(expression_matrix.models_processed, amount_of_models)

    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_r_files,
    )
    @patch(
        "expression_files.expression_files.GeneExpression.is_raw_count_file",
        return_value=False,
    )
    def test_transform_fn_r_format(self, mock_load, mock_is_raw_count_file):
        """
        Assures transform creates data models for r formatted files correctly.
        """
        expression_matrix = DenseIngestor(
            "../tests/data/r_format_text.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        expression_matrix.transform()

    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_genes_batched,
    )
    @patch(
        "expression_files.expression_files.GeneExpression.is_raw_count_file",
        return_value=False,
    )
    def test_transform_fn_batch(self, mock_load, mock_is_raw_count_file):
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
            expression_matrix.transform()

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
            expression_matrix.transform()

    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_no_exp_data,
    )
    @patch(
        "expression_files.expression_files.GeneExpression.is_raw_count_file",
        return_value=False,
    )
    def test_transform_fn_no_exp_data(self, mock_load, mock_is_raw_count_file):
        """
        Confirms gene models are created even when there is no expression data
        """

        expression_matrix = DenseIngestor(
            "../tests/data/dense_matrix_no_exp_data.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
        )
        with patch(
            "expression_files.expression_files.GeneExpression.DATA_ARRAY_BATCH_SIZE",
            new_callable=PropertyMock,
            return_value=4,
        ):
            expression_matrix.transform()
