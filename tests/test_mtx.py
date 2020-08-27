"""Test test_mtx.py

These tests verify:
    - execute_ingest calls proper functions
    - self.genes and self.cells are set properly in
        extract_feature_barcode_matrices()
"""
import unittest
import sys
from unittest.mock import patch, MagicMock, PropertyMock
from bson.objectid import ObjectId

from test_expression_files import mock_load_genes_batched
from mock_data.expression.matrix_mtx.gene_models import mtx_gene_models
from mock_data.expression.matrix_mtx.data_arrays import mtx_data_arrays

sys.path.append("../ingest")
from expression_files.mtx import MTXIngestor
from expression_files.expression_files import GeneExpression

from ingest_files import DataArray


def mock_load_mtx(documents, collection_name):
    """Enables overwriting of GeneExpression.load() with this placeholder.
    GeneExpression.load() is called multiple times. This method will verify
    models in the arguments have the expected values.
    """
    # _id and linear_data_id are unique identifiers and can not be predicted
    # so we exclude it from the comparison
    for document in documents:
        model_name = document["name"]
        if collection_name == GeneExpression.COLLECTION_NAME:
            del document["_id"]
            assert document == mtx_gene_models[model_name]
        if collection_name == DataArray.COLLECTION_NAME:
            del document["linear_data_id"]
            assert document == mtx_data_arrays[model_name]


class TestMTXIngestor(unittest.TestCase):
    def test_check_bundle(self):
        mtx_dimensions = [5, 4, 50]
        expected_genes = mtx_dimensions[0]
        expected_barcodes = mtx_dimensions[1]
        genes = [1, 2, 3, 4, 5]
        barcodes = ["cell0", "cell1", "cell3", "cell4"]
        self.assertTrue(MTXIngestor.check_bundle(barcodes, genes, mtx_dimensions))

        gene_short = [1, 2, 3, 4]
        with self.assertRaises(ValueError) as cm:
            MTXIngestor.check_bundle(barcodes, gene_short, mtx_dimensions)
        expected_msg = (
            f"Expected {expected_barcodes} cells and {expected_genes} genes. "
            f"Got {len(barcodes)} cells and {len(gene_short)} genes."
        )
        self.assertEqual(str(cm.exception), expected_msg)

        barcodes_bad = ["cell0", "cell1", "cell3"]
        with self.assertRaises(ValueError) as cm:
            MTXIngestor.check_bundle(barcodes_bad, genes, mtx_dimensions)
        expected_msg = (
            f"Expected {expected_barcodes} cells and {expected_genes} genes. "
            f"Got {len(barcodes_bad)} cells and {len(genes)} genes."
        )
        self.assertEqual(str(cm.exception), expected_msg)

        with self.assertRaises(ValueError) as cm:
            MTXIngestor.check_bundle(barcodes_bad, gene_short, mtx_dimensions)
        expected_msg = (
            f"Expected {expected_barcodes} cells and {expected_genes} genes. "
            f"Got {len(barcodes_bad)} cells and {len(gene_short)} genes."
        )
        self.assertEqual(str(cm.exception), expected_msg)

    def test_is_sorted(self):
        visited_nums = [0]
        sorted_nums = [2, 2, 2, 3, 4, 4, 4, 5]
        for num in sorted_nums:
            self.assertTrue(MTXIngestor.is_sorted(num, visited_nums))
            if num not in visited_nums:
                visited_nums.append(num)

        visited_nums = [0]
        unsorted_nums = [1, 2, 2, 4, 6, 5]
        truth_values = []
        for num in unsorted_nums:
            truth_values.append(MTXIngestor.is_sorted(num, visited_nums))
            if num not in visited_nums:
                visited_nums.append(num)
        self.assertFalse(all(truth_values))

    def test_get_mtx_dimensions(self):
        file_handler = open("data/AB_toy_data_toy.matrix.mtx")
        dimensions = MTXIngestor.get_mtx_dimensions(file_handler)
        self.assertEqual([80, 272, 4352], dimensions)

        # Dimension contains a letter
        file_handler = open("data/bad_dimensions.mtx.txt")
        self.assertRaises(Exception, MTXIngestor.get_mtx_dimensions, file_handler)

        with self.assertRaises(ValueError) as cm:
            file_handler = open("data/no_data.mtx.txt")
            MTXIngestor.get_mtx_dimensions(file_handler)
        self.assertEqual("MTX file did not contain data", str(cm.exception))

    def test_check_duplicates(self):

        values = ["2", "4", "5", "7"]
        self.assertTrue(MTXIngestor.check_duplicates(values, "scores"))
        dup_values = ["foo1", "f002", "foo1", "foo3"]
        self.assertRaises(
            ValueError, MTXIngestor.check_duplicates, dup_values, "scores"
        )

    @patch("expression_files.expression_files.GeneExpression.check_unique_cells")
    def test_check_valid(self, mock_check_unique_cells):
        """Confirms the errors are correctly promulgated"""
        query_params = (ObjectId("5dd5ae25421aa910a723a337"), MagicMock())

        mock_check_unique_cells.return_value = True
        self.assertTrue(
            MTXIngestor.check_valid(
                ["CELL01", "CELL02", "Cell03"],
                ["gene1", "0", "1"],
                [3, 3, 25],
                query_params,
            )
        )
        duplicate_barcodes = ["CELL01", "CELL02", "CELL02"]
        with self.assertRaises(ValueError) as cm:
            MTXIngestor.check_valid(
                duplicate_barcodes, ["gene1", "0", "1"], [3, 3, 25], query_params
            )

        expected_dup_msg = (
            "Duplicate values are not allowed. "
            "There are 1 duplicates in the barcode file"
        )
        self.assertEqual(expected_dup_msg, str(cm.exception))

        short_genes = ["gene1", "0"]
        # Should concatenate the two value errors
        mock_check_unique_cells.return_value = True
        with self.assertRaises(ValueError) as cm:
            MTXIngestor.check_valid(
                duplicate_barcodes, short_genes, [3, 3, 25], query_params
            )
        expected_msg = (
            "Expected 3 cells and 3 genes. Got 3 cells and 2 genes.;"
            f" {expected_dup_msg}"
        )
        self.assertEqual(expected_msg, str(cm.exception))

    @patch("expression_files.expression_files.GeneExpression.load")
    @patch(
        "expression_files.mtx.MTXIngestor.transform",
        return_value=[({"foo1": "foo2"}, "name")],
    )
    @patch(
        "expression_files.expression_files.GeneExpression.check_unique_cells",
        return_value=True,
    )
    def test_execute_ingest(self, mock_load, mock_transform, mock_has_unique_cells):
        """
            Integration test for execute_ingest()
        """

        expression_matrix = MTXIngestor(
            "../tests/data/matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/AB_toy_data_toy.barcodes.tsv",
        )
        # When is_valid_format() is false exception should be raised
        with self.assertRaises(ValueError) as error:
            expression_matrix.execute_ingest()
        self.assertEqual(
            str(error.exception),
            "Expected 25 cells and 33694 genes. Got 272 cells and 80 genes.",
        )
        expression_matrix = MTXIngestor(
            "../tests/data/AB_toy_data_toy.matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/AB_toy_data_toy.barcodes.tsv",
        )

        expression_matrix.execute_ingest()
        self.assertTrue(mock_transform.called)

    def test_unsorted_mtx_transform(self):
        """
        Tests if value error is raised when mtx file is unsorted,
        """
        expression_matrix = MTXIngestor(
            "../tests/data/unsorted_mtx.mtx.txt",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/AB_toy_data_toy.barcodes.tsv",
        )
        expression_matrix.extract_feature_barcode_matrices()
        with self.assertRaises(ValueError) as cm:
            expression_matrix.transform()
        self.assertEqual("MTX file must be sorted", str(cm.exception))

    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_mtx,
    )
    def test_transform_fn(self, mock_load):
        """
        Assures transform function creates gene data model correctly
        """
        expression_matrix = MTXIngestor(
            "../tests/data/AB_toy_data_toy.matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/AB_toy_data_toy.barcodes.tsv",
        )
        expression_matrix.extract_feature_barcode_matrices()
        expression_matrix.transform()

    @patch(
        "expression_files.expression_files.GeneExpression.load",
        side_effect=mock_load_genes_batched,
    )
    def test_transform_fn_batch(self, mock_load):
        """
        Assures transform function batches data array creation correctly and
        """
        GeneExpression.DATA_ARRAY_BATCH_SIZE = 7
        expression_matrix = MTXIngestor(
            "../tests/data/AB_toy_data_toy.matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/AB_toy_data_toy.barcodes.tsv",
        )
        with patch(
            "expression_files.expression_files.GeneExpression.DATA_ARRAY_BATCH_SIZE",
            new_callable=PropertyMock,
            return_value=7,
        ):
            expression_matrix.extract_feature_barcode_matrices()
            expression_matrix.transform()
