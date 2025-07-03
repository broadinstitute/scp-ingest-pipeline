
import unittest
import sys
import os
import json
import gzip
import pandas as pd
from bson.objectid import ObjectId

sys.path.append("../ingest")
from dot_plot_genes import DotPlotGenes


class TestDotPlotGenes(unittest.TestCase):
    TEST_PREFIX = 'dot_plot_genes_test'
    study_id = ObjectId()
    study_file_id = ObjectId()

    def setup_dense(self):
        return DotPlotGenes(
            study_id=self.study_id,
            study_file_id=self.study_file_id,
            cluster_group_id=f"{self.TEST_PREFIX}_{ObjectId()}",
            annotation_file='data/metadata_example.txt',
            matrix_file_path='data/dense_expression_matrix.txt',
            matrix_file_type='dense',
            cluster_file='data/cluster_example.txt',
        )

    def setup_sparse(self):
        return DotPlotGenes(
            study_id=self.study_id,
            study_file_id=self.study_file_id,
            cluster_group_id=f"{self.TEST_PREFIX}_{ObjectId()}",
            annotation_file='data/mtx/metadata_mtx_barcodes.tsv',
            matrix_file_path='data/mtx/matrix_with_header.mtx',
            gene_file='data/mtx/sampled_genes.tsv',
            barcode_file='data/mtx/barcodes.tsv',
            matrix_file_type='mtx',
            cluster_file='data/mtx/cluster_mtx_barcodes.tsv'
        )

    def test_set_annotation_map(self):
        dot_plot = self.setup_dense()
        dot_plot.set_annotation_map()
        self.assertEqual(['Cluster--group--study', 'Sub-Cluster--group--study'], list(dot_plot.annotation_map.keys()))
        expected_cells = set(['CELL_0001', 'CELL_0002', 'CELL_0003', 'CELL_0004', 'CELL_0005'])
        cluster_cells = [f"CELL_000{i}" for i in range(1, 16)]
        self.assertEqual(expected_cells, dot_plot.annotation_map['Cluster--group--study']['CLST_A'])
        self.assertEqual(cluster_cells, dot_plot.cluster_cells)

    def test_render_expression_dense(self):
        dot_plot = self.setup_dense()
        dot_plot.render_gene_expression()
        expected_data = json.loads(open(f"data/expression_writer/gene_dicts/Sergef.json").read())
        rendered_data = json.loads(gzip.open(f"{dot_plot.cluster_name}/Sergef.json").read())
        self.assertEqual(expected_data, rendered_data)