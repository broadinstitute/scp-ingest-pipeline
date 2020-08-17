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
    def test_extract_feature_barcode_matrices(self):
        self.ingestor.extract_feature_barcode_matrices()

        amount_of_genes = int(self.ingestor.mtx_description[0])
        amount_of_cells = int(self.ingestor.mtx_description[1])
        self.assertEqual(len(self.ingestor.genes), amount_of_genes)
        self.assertEqual(len(self.ingestor.cells), amount_of_cells)

    def test_get_matrix_dimensions(self):
        file_handler = open("data/AB_toy_data_toy.matrix.mtx")
        deminsions = MTXIngestor.get_mtx_dimensions(file_handler)
        self.assertEqual([80, 272, 4352], deminsions)

    def test_check_duplicate_genes(self):

    def test_is_sorted(self):
        visited_nums = [0]
        sorted_nums = [1, 1, 2, 3, 4, 4, 4, 5]
        for num in sorted_nums:
            self.assertTrue(MTXIngestor.is_sorted(num, visited_nums))
            if num not in visited_nums:
                visited_nums.append(num)

        visited_nums = [0]
        unsorted_nums = [1, 2, 2, 3, 8, 9]
        truth_values = []
        self.assertRaises(ValueError, MTXIngestor.is_sorted, num, visited_nums)

    # Make sure lines are length of mtx_descriptor
    # @patch(MTXIngestor.transform)
    # @patch('mtx.MTXIngestor.extract_feature_barcode_matrices')
    # @patch('mtx.close')
    # def test_execute_ingest(self, mock_transform):
    #     # mock_extract_feature_barcode_matrices):
    #     ingestor = MTXIngestor('data/AB_toy_data_toy.matrix.mtx',
    #     '5dd5ae25421aa910a723a337',
    #     '5d276a50421aa9117c982845',
    #     gene_file= 'data/AB_toy_data_toy.genes.tsv',
    #     barcode_file='data/AB_toy_data_toy.barcodes.tsv')
    #     print(ingestor.mtx_desciption)
    #     ingestor.execute_ingest()
    #     mock_transform.assert_called()
    # mock_extract_feature_barcode_matrices.assert_called()
