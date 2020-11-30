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


sys.path.append("../ingest")
from expression_files.mtx import MTXIngestor
from expression_files.expression_files import GeneExpression
from test_expression_files import mock_load_genes_batched, mock_expression_load
from ingest_files import IngestFiles


class TestMTXIngestor(unittest.TestCase):
    GeneExpression.load = mock_expression_load

    def test_get_features(self):
        single_column = "0610007P14Rik"
        self.assertEqual(
            MTXIngestor.get_features(single_column), ("0610007P14Rik", "0610007P14Rik")
        )
        two_columns = "ENSMUSG00000051951\tXkr4"
        self.assertEqual(
            MTXIngestor.get_features(two_columns), ("ENSMUSG00000051951", "Xkr4")
        )

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

    def test_sort_mtx(self):
        import filecmp
        import os

        expected_sorted_mtx = "../tests/data/mtx/sorted_matrix.mtx"
        unsorted_mtx = "../tests/data/mtx/unsorted_matrix.mtx"
        zipped_unsorted_mtx = "../tests/data/mtx/unsorted_matrix.mtx.gz"

        sorted_mtx = MTXIngestor.sort_mtx(unsorted_mtx, open(unsorted_mtx))

        # Verify files have the same contents
        self.assertTrue(filecmp.cmp(sorted_mtx, expected_sorted_mtx))
        # Delete sorted MTX file
        os.remove(sorted_mtx)

        zipped_ingest_file = IngestFiles(
            zipped_unsorted_mtx, ["text/tab-separated-values"]
        )
        zipped_file = zipped_ingest_file.resolve_path(zipped_unsorted_mtx)[0]
        sorted_mtx = MTXIngestor.sort_mtx(zipped_unsorted_mtx, zipped_file)
        # Verify files have the same contents
        self.assertTrue(filecmp.cmp(sorted_mtx, expected_sorted_mtx))
        os.remove(sorted_mtx)

    def test_get_data_start_line_number(self):
        mtx_file_handler = open("data/mtx/unsorted_mtx.mtx.txt")
        self.assertEqual(3, MTXIngestor.get_data_start_line_number(mtx_file_handler))

        mtx_file_handler = open("data/mtx/bad_format_has_character.mtx")
        self.assertRaises(
            ValueError, MTXIngestor.get_data_start_line_number, mtx_file_handler
        )

        mtx_file_handler = open("data/mtx/bad_format_has_space.mtx")
        self.assertRaises(
            IndexError, MTXIngestor.get_data_start_line_number, mtx_file_handler
        )

        # Test for empty file
        empty_file_handler = open("data/empty_file.txt")
        self.assertRaises(
            ValueError, MTXIngestor.get_data_start_line_number, empty_file_handler
        )

    def test_is_sorted(self):
        self.assertTrue(
            MTXIngestor.is_sorted(
                "data/mtx/AB_toy_data_toy.matrix.mtx",
                open("data/mtx/AB_toy_data_toy.matrix.mtx"),
            )
        )
        self.assertFalse(
            MTXIngestor.is_sorted(
                "data/mtx/unsorted_mtx.mtx.txt", open("data/mtx/unsorted_mtx.mtx.txt")
            )
        )

        # Test empty file
        self.assertRaises(
            ValueError,
            MTXIngestor.is_sorted,
            "data/empty_file.txt",
            open("data/empty_file.txt"),
        )

    def test_get_mtx_dimensions(self):
        file_handler = open("data/mtx/AB_toy_data_toy.matrix.mtx")
        dimensions = MTXIngestor.get_mtx_dimensions(file_handler)
        self.assertEqual([80, 272, 4352], dimensions)

        # Dimension contains a letter
        file_handler = open("data/mtx/bad_dimensions.mtx.txt")
        self.assertRaises(Exception, MTXIngestor.get_mtx_dimensions, file_handler)

        with self.assertRaises(ValueError) as cm:
            file_handler = open("data/mtx/no_data.mtx.txt")
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
    @patch(
        "expression_files.expression_files.GeneExpression.is_raw_count_file",
        return_value=False,
    )
    def test_execute_ingest(
        self, mock_load, mock_transform, mock_has_unique_cells, mock_is_raw_count_file
    ):
        """
            Integration test for execute_ingest()
        """

        expression_matrix = MTXIngestor(
            "../tests/data/mtx/matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        )
        # When is_valid_format() is false exception should be raised
        with self.assertRaises(ValueError) as error:
            expression_matrix.execute_ingest()
        self.assertEqual(
            str(error.exception),
            "Expected 25 cells and 33694 genes. Got 272 cells and 80 genes.",
        )
        expression_matrix = MTXIngestor(
            "../tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        )

        expression_matrix.execute_ingest()
        self.assertTrue(mock_transform.called)

    def test_transform_fn(self):
        """
        Assures transform function creates gene data model correctly
        """
        expression_matrix = MTXIngestor(
            "../tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        )
        expression_matrix.test_models = None
        expression_matrix.models_processed = 0
        expression_matrix.extract_feature_barcode_matrices()
        expression_matrix.is_raw_count = False
        expression_matrix.transform()
        amount_of_models = len(
            expression_matrix.test_models["data_arrays"].keys()
        ) + len(expression_matrix.test_models["gene_models"].keys())
        self.assertEqual(expression_matrix.models_processed, amount_of_models)

        # MTX with a single column feature file
        expression_matrix = MTXIngestor(
            "../tests/data/mtx/one_column_feature_file_data_models.mtx",
            "5f2f58eba0845f8c4cf1dc12",
            "5dd5ae25421aa910a723a447",
            gene_file="../tests/data/mtx/one_column_feature_file_data_models.genes.tsv",
            barcode_file="../tests/data/mtx/one_column_feature_file_data_models.barcodes.tsv",
        )
        expression_matrix.test_models = None
        expression_matrix.models_processed = 0
        expression_matrix.extract_feature_barcode_matrices()
        expression_matrix.is_raw_count = False
        expression_matrix.transform()
        amount_of_models = len(expression_matrix.test_models["data_arrays"]) + len(
            expression_matrix.test_models["gene_models"]
        )
        self.assertEqual(expression_matrix.models_processed, amount_of_models)

        # Checks models for unsorted MTX file
        with patch(
            "expression_files.expression_files.GeneExpression.check_unique_cells",
            return_value=True,
        ), patch(
            "expression_files.expression_files.GeneExpression.is_raw_count_file",
            return_value=False,
        ):
            expression_matrix = MTXIngestor(
                "../tests/data/mtx/AB_toy_data_toy.unsorted_mtx.mtx",
                "5d276a50421aa9117c982845",
                "5dd5ae25421aa910a723a337",
                gene_file="../tests/data/mtx/AB_toy_data_toy.genes.tsv",
                barcode_file="../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
            )
            expression_matrix.test_models = None
            expression_matrix.models_processed = 0
            expression_matrix.execute_ingest()
            amount_of_models = len(expression_matrix.test_models["data_arrays"]) + len(
                expression_matrix.test_models["gene_models"]
            )
            self.assertEqual(expression_matrix.models_processed, amount_of_models)

        # Checks models for raw count matrices
        expression_matrix = MTXIngestor(
            "../tests/data/mtx/raw_count.sparse.mtx",
            "5d276a50421aa9117c982845",
            "5dd5ae25421aa910a723a337",
            gene_file="../tests/data/mtx/raw_count.genes.tsv",
            barcode_file="../tests/data/mtx/raw_count.barcodes.tsv",
        )
        expression_matrix.is_raw_count = True
        expression_matrix.test_models = None
        expression_matrix.models_processed = 0
        expression_matrix.extract_feature_barcode_matrices()
        expression_matrix.transform()
        amount_of_models = len(expression_matrix.test_models["data_arrays"])
        self.assertEqual(expression_matrix.models_processed, amount_of_models)

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
            "../tests/data/mtx/AB_toy_data_toy.matrix.mtx",
            "5dd5ae25421aa910a723a447",
            "5f2f58eba0845f8c4cf1dc12",
            gene_file="../tests/data/mtx/AB_toy_data_toy.genes.tsv",
            barcode_file="../tests/data/mtx/AB_toy_data_toy.barcodes.tsv",
        )
        with patch(
            "expression_files.expression_files.GeneExpression.DATA_ARRAY_BATCH_SIZE",
            new_callable=PropertyMock,
            return_value=7,
        ):
            expression_matrix.extract_feature_barcode_matrices()
            expression_matrix.transform()
