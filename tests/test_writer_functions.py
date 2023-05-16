import sys
import unittest
import uuid
import os
import json
import gzip

sys.path.append('../ingest')
from writer_functions import (
    is_gz_file,
    round_exp,
    encode_cluster_name,
    open_file,
    make_data_dir,
    get_matrix_size,
    get_cluster_cells,
    load_entities_as_list,
    process_sparse_fragment,
    extract_sparse_line,
    process_dense_line,
    filter_expression_for_cluster,
    write_gene_scores,
)


class TestWriterFunctions(unittest.TestCase):
    def test_is_gz_file(self):
        self.assertTrue(is_gz_file('data/mtx/unsorted_matrix.mtx.gz'))
        self.assertFalse(is_gz_file('data/mtx/matrix.mtx'))

    def test_round_exp(self):
        self.assertEqual(1.123, round_exp(1.123456789, 3))

    def test_encode_cluster_name(self):
        name = 'Cluster with spaces'
        expected_name = 'Cluster_with_spaces'
        self.assertEqual(expected_name, encode_cluster_name(name))
        name_with_pos = 'CD4+ T cells'
        expected_name_pos = 'CD4pos_T_cells'
        self.assertEqual(expected_name_pos, encode_cluster_name(name_with_pos))

    def test_open_file(self):
        filepath = 'data/mtx/unsorted_matrix.mtx.gz'
        file_io, pathname = open_file(filepath)
        self.assertEqual(filepath, pathname)
        # assert correct IO reader was returned by reading line
        expected_line = "%%MatrixMarket matrix coordinate integer general\n"
        self.assertEqual(expected_line, file_io.readline())

    def test_make_data_dir(self):
        dir_name = f"{uuid.uuid4()}"
        make_data_dir(dir_name)
        self.assertTrue(os.path.exists(dir_name))
        self.assertTrue(os.path.exists(f"{dir_name}/gene_entries"))

    def test_get_matrix_size(self):
        plain_mtx = 'data/mtx/unsorted_matrix.mtx'
        self.assertEqual(165, get_matrix_size(plain_mtx))
        # checking the size of gzipped files will give the _uncompressed_ size
        gzipped_mtx = 'data/mtx/unsorted_matrix.mtx.gz'
        self.assertEqual(165, get_matrix_size(gzipped_mtx))

    def test_get_cluster_cells(self):
        expected_cells = list(f"CELL_000{i}" for i in range(1, 16))
        self.assertEqual(expected_cells, get_cluster_cells('data/cluster_example.txt'))

    def test_load_entities_as_list(self):
        with open('data/mtx/barcodes.tsv') as file:
            barcodes = load_entities_as_list(file)
            self.assertEqual(25, len(barcodes))
            self.assertIn('GACTACAGTAACGCGA-1', barcodes)

    def test_process_sparse_fragment(self):
        barcode_file = open('data/mtx/barcodes.tsv')
        barcodes = list(map(str.strip, barcode_file.readlines()))
        data_dir = 'data/writer_functions'
        # barcodes & cluster cells should be identical in this example
        process_sparse_fragment('OXCT2__entries.txt', barcodes, barcodes, data_dir)
        self.assertTrue(os.path.exists(f"{data_dir}/OXCT2.json"))
        rendered_data = json.loads(gzip.open(f"{data_dir}/OXCT2.json").read())
        expected_data = json.loads(open(f"{data_dir}/OXCT2.orig.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_extract_sparse_line(self):
        line = "1 2 3.456789\n"
        expected_data = [1, 2, 3.456]
        self.assertTrue(expected_data, extract_sparse_line(line))

    def test_process_dense_line(self):
        line = "\"Gad1\"\t0\t0\t1.23423\t0\t0\n"
        cluster_cells = [
            'CTAACTTGTTCCATGA-1',
            'CTCCTAGGTCTCATCC-1',
            'CTCGGAGTCGTAGGAG-1',
            'CTGAAACAGGGAAACA-1',
            'GACTACAGTAACGCGA-1',
        ]
        matrix_cells = [
            'CTCGGAGTCGTAGGAG-1',
            'CTGAAACAGGGAAACA-1',
            'CTCCTAGGTCTCATCC-1',
            'GACTACAGTAACGCGA-1',
            'CTAACTTGTTCCATGA-1',
        ]
        expected_data = [0, 1.234, 0, 0, 0]
        data_dir = 'data/writer_functions'
        process_dense_line(line, matrix_cells, cluster_cells, data_dir)
        self.assertTrue(os.path.exists(f"{data_dir}/Gad1.json"))
        rendered_data = json.loads(gzip.open(f"{data_dir}/Gad1.json").read())
        self.assertEqual(expected_data, rendered_data)

    def test_filter_expression_for_cluster(self):
        cluster_cells = [
            'CTAACTTGTTCCATGA-1',
            'CTCCTAGGTCTCATCC-1',
            'CTCGGAGTCGTAGGAG-1',
            'CTGAAACAGGGAAACA-1',
            'GACTACAGTAACGCGA-1',
        ]
        matrix_cells = [
            'GACTACAGTAACGCGA-1',
            'CTGAAACAGGGAAACA-1',
            'GCATACAGTACCGTTA-1',
            'GCGAGAACAAGAGGCT-1',
            'GCTTCCACAGCTGTAT-1',
        ]
        expression = [1.12345, 2.23456, 3.34567, 4.45678, 5.56789]
        expected_data = [0, 0, 0, 2.23456, 1.12345]
        # filter_expression_for_cluster returns a generator, so cast as list
        filtered_data = list(
            filter_expression_for_cluster(cluster_cells, matrix_cells, expression)
        )
        self.assertEqual(expected_data, filtered_data)

    def test_write_gene_scores(self):
        data = [1.234, 2.345, 3.456, 4.567, 5.678, 6.789]
        gene = 'Egfr'
        data_dir = 'data/writer_functions'
        write_gene_scores(gene, data, data_dir)
        rendered_data = json.loads(gzip.open(f"{data_dir}/{gene}.json").read())
        self.assertEqual(data, rendered_data)
