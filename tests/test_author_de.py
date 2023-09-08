"""Test author_de.py

# Run all tests in a manner that shows report_issues output
pytest test_author_de.py -s
"""

import unittest
import sys
import glob
import os
import pytest

sys.path.append("../ingest")
from author_de import AuthorDifferentialExpression


class TestDifferentialExpression(unittest.TestCase):

    def test_execute(self):
        header_refmap = {
            'genes': 'genes',
            'group': 'group',
            'comparison_group': 'comparison_group',
            'size': 'logfoldchanges',
            'significance': 'qval'
        }

        author_de = AuthorDifferentialExpression(
            'cluster_umap_txt',
            'General_Celltype',
            'SCPdev',
            'study',
            'wilcoxon',
            '../tests/data/author_de/lfc_qval_scanpy-like.csv',
            header_refmap
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

    def test_column_ordering(self):
        """Tests that processed output has specified size metric as 1st, significance as 2nd.

        This order of "genes", "<size_metric>", "<significance_metric>" is needed by the UI.
        """

        header_refmap = {
            'genes': 'genes',
            'group': 'group',
            'comparison_group': 'comparison_group',
            'size': 'avg_log2FC',
            'significance': 'p_val_adj'
        }

        author_de = AuthorDifferentialExpression(
            'cluster_umap_txt',
            'General_Celltype',
            'SCPdev',
            'study',
            'wilcoxon',
            '../tests/data/author_de/pval_lfc_pvaladj_seurat-like.tsv',
            header_refmap
        )
        author_de.execute()

        # Test transformation of one-vs-rest author DE
        one_vs_rest_de_output_file = "cluster_umap_txt--General_Celltype--B_cells--study--wilcoxon.tsv"
        with open(one_vs_rest_de_output_file) as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 41, f"Expected 41 files in: {one_vs_rest_de_output_file}")
        expected_line_0 = 'genes	avg_log2FC	p_val_adj	pct.2	pct.1	p_val	cluster'
        self.assertEqual(lines[0].strip(), expected_line_0)
        expected_line_1 = '0	ACE2	0.04469246812	0.561922816	0.2722805917	0.7277194083	0.561922816	6'
        self.assertEqual(lines[1].strip(), expected_line_1)

    def test_size_and_significance_validation(self):
        """Tests validation that file has specified size and significance headers
        """
        with pytest.raises(ValueError) as exc_info:
            header_refmap = {
                'genes': 'genes',
                'group': 'group',
                'comparison_group': 'comparison_group',
                'size': 'OTHER_SIZE',
                'significance': 'OTHER_SIGNIFICANCE'
            }
            author_de = AuthorDifferentialExpression(
                'cluster_umap_txt',
                'General_Celltype',
                'SCPdev',
                'study',
                'wilcoxon',
                '../tests/data/author_de/pval_lfc_pvaladj_seurat-like.tsv',
                header_refmap
            )
            author_de.execute()

        expected_msg = "Column headers must include \"OTHER_SIZE\" and \"OTHER_SIGNIFICANCE\".  No such size or significance metrics found in headers: ['pct.2', 'pct.1', 'p_val_adj', 'p_val', 'cluster', 'avg_log2FC']"
        self.assertEqual(str(exc_info.value), expected_msg)

    def test_seurat_findallmarkers(self):
        """Tests handling for Seurat FindAllMarkers() format
        """

        # The structure (but no specific metric values or genes) of this file
        # matches the first DE file uploaded to SCP by a real user.
        input_filename = 'seurat_findallmarkers_one-vs-rest.csv'

        header_refmap = {
            'genes': 'gene',
            'group': 'cluster',
            'comparison_group': 'None',
            'size': 'avg_log2FC',
            'significance': 'p_val_adj'
        }

        author_de = AuthorDifferentialExpression(
            'cluster_umap_txt',
            'General_Celltype',
            'SCPdev',
            'study',
            'wilcoxon',
            f'../tests/data/author_de/{input_filename}',
            header_refmap
        )
        author_de.execute()

        output_file = "cluster_umap_txt--General_Celltype--10--study--wilcoxon.tsv"
        with open(output_file) as f:
            lines = [line.strip() for line in f.readlines()]

        expected_line_0 = 'genes	avg_log2FC	p_val_adj	pct.2	pct.1	p_val'
        expected_line_1 = '0	CRP	1.082893128	0.0	0.1364	0.3268	0.0'
        expected_line_2 = '1	LGALS3	1.399143463	0.0	0.29848	0.4284	0.0'
        expected_line_3 = '2	PIK3CA	0.7331456728	7.070000000000001e-273	0.34236	0.67932	3.49e-277'

        self.assertEqual(len(lines), 4, f"Expected 4 files in: {output_file}")
        self.assertEqual(lines[0], expected_line_0)
        self.assertEqual(lines[1], expected_line_1)
        self.assertEqual(lines[2], expected_line_2)
        self.assertEqual(lines[3], expected_line_3)

    def teardown_method(self, test_method):
        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        for file in files:
            os.remove(file)



