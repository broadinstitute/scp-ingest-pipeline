"""Test test_mtx.py

These tests verify:
    - execute_ingest calls proper functions
    - self.genes and self.cells are set properly in
        extract_feature_barcode_matrices()
"""
import unittest
import sys


sys.path.append("../ingest")
from expression_files.mtx import MTXIngestor
import sys


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

    def test_get_mtx_dimensions(self):
        file_handler = open("data/AB_toy_data_toy.matrix.mtx")
        dimensions = MTXIngestor.get_mtx_dimensions(file_handler)
        self.assertEqual([80, 272, 4352], dimensions)

    def test_check_duplicate_genes(self):

        values = ["2", "4", "5", "7"]
        self.assertTrue(MTXIngestor.check_duplicate_genes, values)
        dup_values = ["foo1", "f002", "foo1", "foo3"]
        self.assertRaises(ValueError, MTXIngestor.check_duplicate_genes, dup_values)

    def test_is_sorted(self):
        visited_nums = [0]
        sorted_nums = [1, 1, 2, 3, 4, 4, 4, 5]
        for num in sorted_nums:
            self.assertTrue(MTXIngestor.is_sorted(num, visited_nums))
            if num not in visited_nums:
                visited_nums.append(num)

        unsorted = [1, 2, 2, 4]
        self.assertRaises(ValueError, MTXIngestor.is_sorted, unsorted, visited_nums)
