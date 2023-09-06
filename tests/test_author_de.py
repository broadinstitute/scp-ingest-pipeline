"""Test author_de.py

# Run all tests in a manner that shows report_issues output
pytest test_author_de.py -s
"""

import unittest
import sys
import glob
import os

sys.path.append("../ingest")
from author_de import AuthorDifferentialExpression


class TestDifferentialExpression(unittest.TestCase):

    def test_execute(self):
        author_de = AuthorDifferentialExpression(
            'cluster_umap_txt',
            'General_Celltype',
            'SCPdev',
            'study',
            'wilcoxon',
            '../tests/data/author_de/lfc_qval_scanpy-like.csv',
            'logfoldchanges',
            'qval'
        )
        author_de.execute()

        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        self.assertEqual(len(files), 7, 'Expected 7 files')

        # Test transformation of one-vs-rest author DE
        one_vs_rest_de_output_file = "cluster_umap_txt--General_Celltype--B_cells--study--wilcoxon.tsv"
        with open(one_vs_rest_de_output_file) as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 41, f"Expected 41 files in: {one_vs_rest_de_output_file}")
        expected_line_0 = 'genes	logfoldchanges	qval	mean'
        self.assertEqual(lines[0].strip(), expected_line_0)
        expected_line_1 = '0	ACE2	0.9852699170700532	0.146924681159184	0.2277194082712261'

        # Test transformation of pairwise author DE
        pairwise_de_output_file = "cluster_umap_txt--General_Celltype--B_cells--CSN1S1_macrophages--study--wilcoxon.tsv"
        with open(pairwise_de_output_file) as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 41, f"Expected 41 files in: {pairwise_de_output_file}")
        expected_line_0 = 'genes	logfoldchanges	qval	mean'
        self.assertEqual(lines[0].strip(), expected_line_0)
        expected_line_1 = '0	ACE2	0.685269917070053	0.446924681159184	0.727719408271226'
        self.assertEqual(lines[1].strip(), expected_line_1)

    def test_basic_size_and_significance(self):
        """Tests validation that file has "logfoldchanges" and "qval" columns

        TODO: Revamp this test once UI supports flexible columns
        """
        author_de = AuthorDifferentialExpression(
            'cluster_umap_txt',
            'General_Celltype',
            'SCPdev',
            'study',
            'wilcoxon',
            '../tests/data/author_de/pval_lfc_pvaladj_seurat-like.tsv',
            'avg_log2FC',
            'p_val_adj'
        )

        try:
            author_de.execute()
        except ValueError as e:
            expected_msg = "Column headers must include \"logfoldchanges\" and \"qval\".  No size or significance metrics found in headers: ['p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster']"
            self.assertEqual(str(e), expected_msg)

    def teardown_method(self, test_method):
        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        for file in files:
            os.remove(file)



