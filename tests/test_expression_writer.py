from __future__ import annotations
import sys
import unittest
from unittest.mock import patch
import uuid
import os
import json
import gzip
from shutil import rmtree, copy
import glob

sys.path.append('../ingest')
from expression_writer import ExpressionWriter
from writer_functions import load_entities_as_list, make_data_dir, encode_cluster_name


class TestExpressionWriter(unittest.TestCase):
    TEST_PREFIX = 'exp_writer_test_'

    @staticmethod
    def setup_dense_exp_writer(cluster_name) -> ExpressionWriter:
        return ExpressionWriter(
            matrix_file_path='data/dense_expression_matrix.txt',
            matrix_file_type='dense',
            cluster_file_path='data/cluster_example.txt',
            gene_file=None,
            barcode_file=None,
            cluster_name=cluster_name,
        )

    @staticmethod
    def setup_sparse_exp_writer(cluster_name) -> ExpressionWriter:
        return ExpressionWriter(
            matrix_file_path='data/mtx/matrix_with_header.mtx',
            gene_file='data/mtx/sampled_genes.tsv',
            barcode_file='data/mtx/barcodes.tsv',
            matrix_file_type='mtx',
            cluster_file_path='data/mtx/cluster_mtx_barcodes.tsv',
            cluster_name=cluster_name,
        )

    @staticmethod
    def seed_test_gene_entries(data_dir):
        for file in os.listdir('data/expression_writer/gene_entries'):
            copy(
                f"data/expression_writer/gene_entries/{file}",
                f"{data_dir}/gene_entries",
            )

    @classmethod
    def teardown_class(cls):
        logs = glob.glob('expression_scatter_data_*_log.txt')
        for log in logs:
            os.remove(log)
        test_dirs = glob.glob(f"{TestExpressionWriter.TEST_PREFIX}*")
        for dirname in test_dirs:
            rmtree(dirname)

    def test_process_dense_matrix(self):
        cluster_name = encode_cluster_name(f"{self.TEST_PREFIX}dense_{uuid.uuid4()}")
        exp_writer = self.setup_dense_exp_writer(cluster_name)
        exp_writer.render_artifacts()
        self.assertTrue(os.path.exists(cluster_name))
        self.assertTrue(os.path.exists(f"{cluster_name}/Sergef.json"))
        self.assertTrue(os.path.exists(f"{cluster_name}/Itm2a.json"))
        expected_data = json.loads(open(f"data/expression_writer/Sergef.json").read())
        rendered_data = json.loads(gzip.open(f"{cluster_name}/Sergef.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_process_sparse_matrix(self):
        cluster_name = encode_cluster_name(f"{self.TEST_PREFIX}sparse_{uuid.uuid4()}")
        exp_writer = self.setup_sparse_exp_writer(cluster_name)
        exp_writer.render_artifacts()
        self.assertTrue(os.path.exists(cluster_name))
        genes = load_entities_as_list(open(exp_writer.gene_file))
        genes.remove('HOMER2')  # doesn't render gene entry file
        for gene in genes:
            self.assertTrue(os.path.exists(f"{cluster_name}/{gene}.json"))
            self.assertTrue(
                os.path.exists(f"{cluster_name}/gene_entries/{gene}__entries.txt")
            )
        expected_data = json.loads(
            open(f"data/writer_functions/OXCT2.orig.json").read()
        )
        rendered_data = json.loads(gzip.open(f"{cluster_name}/OXCT2.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_get_bucket_name(self):
        with patch('ingest_files.IngestFiles.verify_file_exists', return_value=None):
            with patch(
                'ingest_files.IngestFiles.download_from_bucket',
                return_value='data/mtx/matrix.mtx',
            ):
                bucket_name = f"gs://{uuid.uuid4()}"
                exp_writer = ExpressionWriter(
                    matrix_file_path=f"{bucket_name}/matrix.mtx",
                    cluster_file_path=f"{bucket_name}/cluster.txt",
                    cluster_name='test',
                    gene_file=None,
                    barcode_file=None,
                    matrix_file_type="dense",
                )
                self.assertEqual(bucket_name, exp_writer.get_storage_bucket_name())

    def test_get_file_seek_points(self):
        cluster_name = f"{self.TEST_PREFIX}seek_points_{uuid.uuid4()}"
        exp_writer = self.setup_dense_exp_writer(cluster_name)
        seek_points = exp_writer.get_file_seek_points()
        # note: this is dependent on the number of cores, and depending on your architecture this may differ
        # this test covers cases for both 3 and 4 cores utilized
        if exp_writer.num_cores == 3:
            expected_points = [[161, 264], [265, 265]]
        else:
            expected_points = [[161, 196], [197, 264], [265, 265]]
        self.assertEqual(expected_points, seek_points)

    def test_divide_sparse_matrix(self):
        cluster_name = encode_cluster_name(
            f"{self.TEST_PREFIX}sparse_divide_{uuid.uuid4()}"
        )
        exp_writer = self.setup_sparse_exp_writer(cluster_name)
        make_data_dir(cluster_name)
        genes = load_entities_as_list(open(exp_writer.gene_file))
        exp_writer.divide_sparse_matrix(genes, cluster_name)
        genes.remove('HOMER2')  # doesn't render gene entry file
        for gene in genes:
            self.assertTrue(
                os.path.exists(f"{cluster_name}/gene_entries/{gene}__entries.txt")
            )

    def test_read_sparse_matrix_slice(self):
        cluster_name = encode_cluster_name(
            f"{self.TEST_PREFIX}sparse_slice_{uuid.uuid4()}"
        )
        indexes = [136, 219]
        exp_writer = self.setup_sparse_exp_writer(cluster_name)
        make_data_dir(cluster_name)
        genes = load_entities_as_list(open(exp_writer.gene_file))
        exp_writer.read_sparse_matrix_slice(indexes, genes, cluster_name)
        self.assertTrue(
            os.path.exists(f"{cluster_name}/gene_entries/OXCT2__entries.txt")
        )
        expected_data = open(
            'data/expression_writer/gene_entries/OXCT2__entries.txt'
        ).read()
        rendered_data = open(f"{cluster_name}/gene_entries/OXCT2__entries.txt").read()
        self.assertEqual(expected_data, rendered_data)

    def test_process_sparse_data_fragments(self):
        cluster_name = encode_cluster_name(
            f"{self.TEST_PREFIX}sparse_fragments_{uuid.uuid4()}"
        )
        exp_writer = self.setup_sparse_exp_writer(cluster_name)
        make_data_dir(cluster_name)
        self.seed_test_gene_entries(cluster_name)
        cells = load_entities_as_list(open(exp_writer.barcode_file))
        genes = load_entities_as_list(open(exp_writer.gene_file))
        exp_writer.process_sparse_data_fragments(cells, cells, cluster_name)
        genes.remove('HOMER2')
        for gene in genes:
            self.assertTrue(os.path.exists(f"{cluster_name}/{gene}.json"))

    def test_write_empty_sparse_genes(self):
        cluster_name = encode_cluster_name(
            f"{self.TEST_PREFIX}sparse_empty_genes_{uuid.uuid4()}"
        )
        exp_writer = self.setup_sparse_exp_writer(cluster_name)
        make_data_dir(cluster_name)
        self.seed_test_gene_entries(cluster_name)
        num_cells = 25
        genes = load_entities_as_list(open(exp_writer.gene_file))
        exp_writer.write_empty_sparse_genes(genes, num_cells, cluster_name)
        # only empty gene should be HOMER2
        gene = 'HOMER2.json'
        self.assertTrue(os.path.exists(f"{cluster_name}/{gene}"))
        expected_data = list(0 for _ in range(0, 25))
        empty_data = json.loads(gzip.open(f"{cluster_name}/{gene}").read())
        self.assertEqual(expected_data, empty_data)

    def test_read_dense_matrix_slice(self):
        cluster_name = encode_cluster_name(
            f"{self.TEST_PREFIX}dense_slice_{uuid.uuid4()}"
        )
        exp_writer = self.setup_dense_exp_writer(cluster_name)
        make_data_dir(cluster_name)
        indexes = [161, 264]
        cells = list(f"CELL_000{i}" for i in range(1, 16))
        exp_writer.read_dense_matrix_slice(indexes, cells, cells, cluster_name)
        self.assertTrue(os.path.exists(f"{cluster_name}/Sergef.json"))
        self.assertTrue(os.path.exists(f"{cluster_name}/Itm2a.json"))
