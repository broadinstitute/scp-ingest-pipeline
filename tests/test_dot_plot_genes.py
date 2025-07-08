import unittest
import sys
import os
import json
import gzip
import glob
from shutil import rmtree
from bson.objectid import ObjectId

sys.path.append("../ingest")
from dot_plot_genes import DotPlotGenes


class TestDotPlotGenes(unittest.TestCase):
    TEST_PREFIX = 'dot_plot_genes_test'
    study_id = str(ObjectId())
    study_file_id = str(ObjectId())
    SERGEF_METRICS = {
        "Cluster--group--study": {
            "CLST_A": [4.281, 0.6],
            "CLST_B": [2.869, 0.4],
            "CLST_C": [3.876, 0.6]
        },
        "Sub-Cluster--group--study": {
            "CLST_A_1": [2.364, 0.3333],
            "CLST_A_2": [7.157, 1.0],
            "CLST_B_1": [0.0, 0.0],
            "CLST_B_2": [4.782, 0.6667],
            "CLST_C_1": [4.59, 0.6667],
            "CLST_C_2": [2.805, 0.5]
        }
    }
    OXCT2_METRICS = {
        "Cluster--group--study": {
            "CLST_A": [0.0, 0.0],
            "CLST_B": [0.0, 0.0],
            "CLST_C": [0.525, 0.1333]
        },
        "Sub-Cluster--group--study": {
            "CLST_A_1": [0.0, 0.0],
            "CLST_A_2": [0.0, 0.0],
            "CLST_B_1": [0.0, 0.0],
            "CLST_B_2": [0.0, 0.0],
            "CLST_C_1": [0.274, 0.1667],
            "CLST_C_2": [0.0, 0.0],
            "CLST_C_3": [0.0, 0.0],
            "CLST_C_4": [3.118, 0.5]
        }
    }

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

    @classmethod
    def teardown_class(cls):
        logs = glob.glob('expression_scatter_data_*_log.txt')
        for log in logs:
            os.remove(log)
        test_dirs = glob.glob(f"cluster_{TestDotPlotGenes.TEST_PREFIX}*")
        for dirname in test_dirs:
            rmtree(dirname)

    def dense_cluster_cells(self):
        return [f"CELL_000{i}" for i in range(1, 16)]

    def test_set_annotation_map(self):
        dot_plot = self.setup_dense()
        dot_plot.set_annotation_map()
        self.assertEqual(['Cluster--group--study', 'Sub-Cluster--group--study'], list(dot_plot.annotation_map.keys()))
        expected_cells = set(['CELL_0001', 'CELL_0002', 'CELL_0003', 'CELL_0004', 'CELL_0005'])
        self.assertEqual(expected_cells, dot_plot.annotation_map['Cluster--group--study']['CLST_A'])
        self.assertEqual(self.dense_cluster_cells(), dot_plot.cluster_cells)

    def test_render_expression_dense(self):
        dot_plot = self.setup_dense()
        dot_plot.render_gene_expression()
        expected_data = json.loads(open("data/expression_writer/gene_dicts/Sergef.json").read())
        rendered_data = json.loads(gzip.open(f"{dot_plot.cluster_name}/Sergef.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_render_expression_sparse(self):
        dot_plot = self.setup_sparse()
        dot_plot.render_gene_expression()
        expected_data = json.loads(open("data/expression_writer/gene_dicts/OXCT2.json").read())
        rendered_data = json.loads(gzip.open(f"{dot_plot.cluster_name}/OXCT2.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_compute_pct_exp(self):
        gene_doc = json.loads(open("data/expression_writer/gene_dicts/Sergef.json").read())
        self.assertEqual(0.5333, DotPlotGenes.pct_expression(gene_doc, self.dense_cluster_cells()))

    def test_compute_scaled_mean(self):
        gene_doc = json.loads(open("data/expression_writer/gene_dicts/Sergef.json").read())
        self.assertEqual(3.675, DotPlotGenes.scaled_mean_expression(gene_doc, 0.5333))

    def test_get_gene_name(self):
        cluster = "cluster-foo"
        self.assertEqual("Sergef", DotPlotGenes.get_gene_name(f"{cluster}/Sergef.json"))

    def test_get_gene_dict(self):
        dot_plot = self.setup_dense()
        dot_plot.render_gene_expression()
        expected_data = json.loads(open("data/expression_writer/gene_dicts/Sergef.json").read())
        rendered_data = dot_plot.get_gene_dict(f"{dot_plot.cluster_name}/Sergef.json")
        self.assertEqual(expected_data, rendered_data)

    def test_get_filtered_exp(self):
        filtered_cells = ['CELL_0002', 'CELL_0004']
        gene_doc = json.loads(open("data/expression_writer/gene_dicts/Sergef.json").read())
        filtered_exp = DotPlotGenes.filter_expression_for_label(gene_doc, filtered_cells)
        self.assertEqual({'CELL_0002': 7.092, 'CELL_0004': 7.511}, filtered_exp)

    def test_get_expression_metrics(self):
        dot_plot = self.setup_dense()
        dot_plot.preprocess()
        output_path = f"{dot_plot.cluster_name}/dot_plot_genes"
        os.mkdir(output_path)
        gene_doc = json.loads(gzip.open(f"{dot_plot.cluster_name}/Sergef.json").read())
        rendered_exp = DotPlotGenes.get_expression_metrics(gene_doc, dot_plot.annotation_map)
        self.assertEqual(self.SERGEF_METRICS, rendered_exp)

    def test_process_gene(self):
        dot_plot = self.setup_dense()
        dot_plot.preprocess()
        blank_dot_plot_gene = {
            "study_id": dot_plot.study_id,
            "study_file_id": dot_plot.study_file_id,
            "cluster_group_id": dot_plot.cluster_group_id,
            "exp_scores": {}
        }
        output_path = f"{dot_plot.cluster_name}/dot_plot_genes"
        os.mkdir(output_path)
        DotPlotGenes.process_gene(
            f"{dot_plot.cluster_name}/Sergef.json", output_path, blank_dot_plot_gene, dot_plot.annotation_map
        )
        rendered_path = f"{dot_plot.cluster_name}/dot_plot_genes/Sergef.json"
        self.assertTrue(os.path.exists(rendered_path))
        rendered_gene = json.loads(gzip.open(rendered_path).read())
        self.assertEqual(self.SERGEF_METRICS, rendered_gene['exp_scores'])

    def test_process_all_genes(self):
        dot_plot = self.setup_dense()
        dot_plot.preprocess()
        dot_plot.process_all_genes()
        total_genes = os.listdir(f"{dot_plot.cluster_name}/dot_plot_genes")
        self.assertEqual(3, len(total_genes))
        rendered_path = f"{dot_plot.cluster_name}/dot_plot_genes/Sergef.json"
        self.assertTrue(os.path.exists(rendered_path))
        rendered_gene = json.loads(gzip.open(rendered_path).read())
        self.assertEqual(self.SERGEF_METRICS, rendered_gene['exp_scores'])

    def test_run_processor(self):
        dot_plot = self.setup_sparse()
        dot_plot.run_processor()
        total_genes = os.listdir(f"{dot_plot.cluster_name}/dot_plot_genes")
        self.assertEqual(7, len(total_genes))
        rendered_path = f"{dot_plot.cluster_name}/dot_plot_genes/OXCT2.json"
        self.assertTrue(os.path.exists(rendered_path))
        rendered_gene = json.loads(gzip.open(rendered_path).read())
        self.assertEqual(self.OXCT2_METRICS, rendered_gene['exp_scores'])
