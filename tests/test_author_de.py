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
        test_kwargs = {
            'study_accession': 'SCPdev',
            'annotation_name': 'General_Celltype',
            'annotation_type': 'group',
            'annotation_scope': 'study',
            'method': 'wilcoxon',
            'cluster_name': 'cluster_umap_txt',
            'differential_expression_file': '../tests/data/author_de/lfc_qval_scanpy-like.csv'
        }
        author_de = AuthorDifferentialExpression(
            cluster=None,
            cell_metadata=None,
            **test_kwargs,
        )
        author_de.execute()

        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        self.assertEqual(len(files), 7, 'Expected 7 files')

        pairwise_de_file = "cluster_umap_txt--General_Celltype--B_cells--CSN1S1_macrophages--study--wilcoxon.tsv"
        with open(pairwise_de_file) as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 41, f"Expected 41 files in: {pairwise_de_file}")
        expected_line_0 = 'genes	logfoldchanges	qval	mean'
        self.assertEqual(lines[0].strip(), expected_line_0)
        expected_line_1 = '0	ACE2	0.685269917070053	0.685269917070053	0.446924681159184'
        self.assertEqual(lines[1].strip(), expected_line_1)

    @classmethod
    def teardown_class(cls):
        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        for file in files:
            os.remove(file)



