"""expression_writer.py
extract gene-level expression data from dense/sparse matrix files and process in the
context of a cluster file's cells.
this mimics the `expression` data array in visualization API responses

EXAMPLES

dense matrix:
python3 expression_writer.py --matrix-file ../../tests/data/dense_expression_matrix.txt \
                             --cluster-file ../../tests/data/cluster_example.txt \
                             --cluster-name 'Dense Example'

dense matrix with precision override:
python3 expression_writer.py --matrix-file ../../tests/data/dense_expression_matrix.txt \
                             --cluster-file ../../tests/data/cluster_example.txt \
                             --cluster-name 'Dense Example' --precision 1

sparse matrix:
python3 expression_writer.py --matrix-file ../../tests/data/mtx/matrix_with_header.mtx \
                             --genes-file ../../tests/data/mtx/sampled_genes.tsv \
                             --barcodes-file ../../tests/data/mtx/barcodes.tsv \
                             --cluster-file ../../tests/data/mtx/cluster_mtx_barcodes.tsv \
                             --cluster-name 'Sparse Example'
"""

import json
import os
import re
import argparse
import uuid
import gzip
import time
import multiprocessing
from functools import partial

# regex to split on commas and tabs
COMMA_OR_TAB = r"[,\t]"

# regex to split on comma/tab/space
ALL_DELIM = r"[\s,\t]"

# Gzip magic number
GZIP_MAGIC_NUM = b'\x1f\x8b'

class ExpressionWriter:
    precision = 3
    num_cores = multiprocessing.cpu_count() - 1

    def __init__(
        self, matrix_file_path, cluster_file_path, cluster_name, genes_file, barcodes_file, num_digits, cores
    ):
        self.matrix_file_path = matrix_file_path
        self.cluster_file_path = cluster_file_path
        self.cluster_name = cluster_name
        self.genes_file = genes_file
        self.barcodes_file = barcodes_file
        if num_digits is not None:
            self.precision = num_digits
        if cores is not None:
            self.num_cores = cores

    @staticmethod
    def is_gz_file(filepath):
        """
        Determine if a file is gzipped by checking first two bytes against GZIP_MAGIC_NUM

        :param filepath: (String) path to file to check
        :returns (Bool)
        """
        with open(filepath, 'rb') as test_f:
            return test_f.read(2) == GZIP_MAGIC_NUM

    @staticmethod
    def open_file(filepath):
        """
        Open a file with the correct reader

        :param filepath: (String) path to file to open
        :returns: (TextIOWrapper)
        """
        if ExpressionWriter.is_gz_file(filepath):
            return gzip.open(filepath, 'rt')
        else:
            return open(filepath, 'r')

    @staticmethod
    def make_data_dir(name):
        """
        Make a directory to put output/work files in.  Appends a UUID to the end of the directory
        to allow for quick iteration and side-by-side comparison of outputs without manual cleanup

        :param name: (String) name of directory
        :returns: (String)
        """
        dirname = str(f"{name}-{uuid.uuid4()}")
        print(f"creating data directory at {dirname}")
        os.mkdir(dirname)
        os.mkdir(f"{dirname}/gene_entries") # for slicing sparse matrix files
        return dirname

    @staticmethod
    def get_matrix_size(matrix_file_path):
        """
        Get the actual size of a matrix file, whether gzipped or not
        Because gzipped files only use 4 bytes to store the 'original' file size,
        any files larger that 4GB will only store the modulus of 2 ** 32, not the
        actual size of the file.
        Using seek() and tell() will return the last byte position

        :param matrix_file_path: (String) path to matrix file
        :return: (Int)
        """

        if ExpressionWriter.is_gz_file(matrix_file_path):
            with ExpressionWriter.open_file(matrix_file_path) as matrix_file:
                matrix_file.seek(0, 2)
                return matrix_file.tell()
        else:
            return os.stat(matrix_file_path).st_size

    @staticmethod
    def get_cluster_cells(file_path):
        """
        Return a list of cell names from a cluster file (1st column, starting line 3)

        :param file_path: (String) absolute path to cluster_file
        :returns: (List)
        """
        cells = []
        with ExpressionWriter.open_file(file_path) as file:
            file.readline()
            file.readline()
            for row in file:
                cell = re.split(COMMA_OR_TAB, row)[0]
                cells.append(cell)
        return cells

    @staticmethod
    def load_entities_as_list(file, column = None):
        """
        Read an entire 10X feature/barcode file into a list for parsing sparse data

        :param file: (TextIOWrapper) open file object
        :param column: (Integer) specific column to extract from entity file
        :returns: (List)
        """
        print(f"reading entities from {file.name}")
        if not column:
            return list(map(str.rstrip, file.readlines()))
        else:
            return list(re.split(ALL_DELIM, line.strip())[column] for line in file)

    def get_file_seek_points(self):
        """
        Determine start/stop points in a matrix to process in parallel
        Will read in chunks and return a list of start/stop points
        Ensures breaks on newlines

        :return: (List<List>)
        """
        file_size = ExpressionWriter.get_matrix_size(self.matrix_file_path)
        chunk_size = int(file_size / self.num_cores)
        print(f"determining seek points for {self.matrix_file_path} with chunk size {chunk_size}")
        with ExpressionWriter.open_file(self.matrix_file_path) as matrix_file:
            # fast-forward through any comments if this is a sparse matrix
            first_char = matrix_file.read(1)
            header_lines = 3 if first_char == '%' else 1
            for i in range(header_lines):
                matrix_file.readline()
            current_pos = matrix_file.tell()
            current_seek = [current_pos]
            seek_points = []
            for index in range(self.num_cores):
                seek_point = current_pos + chunk_size
                matrix_file.seek(seek_point)
                current_byte = matrix_file.read(1)
                if current_byte == '':  # eof
                    current_seek.append(file_size)
                    seek_points.append(current_seek)
                    break
                while current_byte != "\n":
                    current_byte = matrix_file.read(1)
                    seek_point += 1
                current_seek.append(seek_point)
                seek_points.append(current_seek)
                current_pos = seek_point + 1
                current_seek = [current_pos]
        return seek_points

    def divide_sparse_matrix(self, genes, data_dir):
        """
        Slice a sparse matrix into 1GB chunks and write out individual
        gene-level files in parallel

        :param genes: (List) gene names from features file
        :param data_dir: (String) name out output dir
        """
        print(f"loading sparse data from {self.matrix_file_path}")
        slice_indexes = self.get_file_seek_points()
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(self.read_sparse_matrix_slice,
                            matrix_file_path=self.matrix_file_path, genes=genes, data_dir=data_dir)
        pool.map(processor, slice_indexes)

    @staticmethod
    def read_sparse_matrix_slice(indexes, matrix_file_path, genes, data_dir):
        """
        Read a sparse matrix using start/stop indexes and append to individual gene-level files

        :param indexes: (List) start/stop index points to read from/to
        :param matrix_file_path: (String): path to matrix file
        :param genes: (List) gene names from features file
        :param data_dir: (String) name out output dir
        """
        start_pos, end_pos = indexes
        print(f"reading {matrix_file_path} at index {start_pos}:{end_pos}")
        with ExpressionWriter.open_file(matrix_file_path) as matrix_file:
            current_pos = start_pos
            matrix_file.seek(current_pos)
            while current_pos < end_pos:
                line = matrix_file.readline()
                gene_idx = int(line.split()[0])
                gene_name = genes[gene_idx - 1]
                fragment_path = f"{data_dir}/gene_entries/{gene_name}__entries.txt"
                with open(fragment_path, 'a+') as file:
                    file.write(line)
                current_pos += len(line)

    def process_sparse_data_fragments(self, barcodes, cluster_cells, cluster_name, data_dir):
        """
        Find and process all generated single-gene sparse data fragments

        :param barcodes: (List) list of cell barcodes
        :param cluster_cells: (List) list of cells from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        fragments = os.listdir(f"{data_dir}/gene_entries")
        print(f"subdivision complete, processing {len(fragments)} fragments")
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(self.process_fragment,
                            barcodes=barcodes, cluster_cells=cluster_cells,
                            cluster_name=cluster_name, data_dir=data_dir)
        pool.map(processor, fragments)

    def process_fragment(self, fragment_name, barcodes, cluster_cells, cluster_name, data_dir):
        """
        Process a single-gene sparse matrix fragment and write expression array

        :param fragment_name: (String) name of gene-level fragment file
        :param barcodes: (List) list of cell barcodes
        :param cluster_cells: (List) list of cells from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        gene_name = fragment_name.split('__')[0]
        full_path = f"{data_dir}/gene_entries/{fragment_name}"
        observed_cells = []
        exp_vals = []
        for line in ExpressionWriter.open_file(full_path):
            gene_idx, barcode_idx, exp_val = self.extract_sparse_line(line)
            observed_cells.append(barcodes[barcode_idx - 1])
            exp_vals.append(exp_val)
        filtered_expression = ExpressionWriter.filter_expression_for_cluster(
            cluster_cells, observed_cells, exp_vals
        )
        ExpressionWriter.write_gene_scores(gene_name, cluster_name, filtered_expression, data_dir)

    def extract_sparse_line(self, line):
        """
        Process a single line from a sparse matrix and extract values as integers

        :param line: (String) single line from matrix file
        :return: (List): values as integers
        """
        gene_idx, barcode_idx, raw_exp = line.rstrip().split(' ')
        return [int(gene_idx), int(barcode_idx), round(float(raw_exp), self.precision)]

    def write_empty_sparse_genes(self, genes, num_cluster_cells, cluster_name, data_dir):
        """
        Write out empty arrays of expression values for genes with no significant expression in a sparse matrix

        :param genes: (List) gene names from features file
        :param num_cluster_cells: (Integer) number of cells from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        gene_fragments = filter(lambda file: file[0] != '.', os.listdir(f"{data_dir}/gene_entries"))
        significant_genes = set([gene.split('__')[0] for gene in gene_fragments])
        missing_genes = [gene for gene in genes if gene not in significant_genes]
        empty_expression = [0] * num_cluster_cells
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(ExpressionWriter.write_gene_scores,
                            cluster_name=cluster_name, exp_values=empty_expression, data_dir=data_dir)
        pool.map(processor, missing_genes)

    def process_dense_data(self, cluster_cells, cluster_name, data_dir):
        """
        Main handler to read dense matrix data and process entries at the gene level

        :param cluster_cells: (List) cell names from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        pool = multiprocessing.Pool(self.num_cores)
        slice_indexes = self.get_file_seek_points()
        matrix_file = ExpressionWriter.open_file(self.matrix_file_path)
        header = matrix_file.readline().rstrip()
        values = re.split(COMMA_OR_TAB, header)
        cells = values[1:]
        processor = partial(
            self.read_dense_matrix_slice,
            matrix_cells=cells, cluster_cells=cluster_cells, cluster_name=cluster_name, data_dir=data_dir
        )
        pool.map(processor, slice_indexes)

    def read_dense_matrix_slice(self, indexes, matrix_cells, cluster_cells, cluster_name, data_dir):
        """
        Read a dense matrix using start/stop indexes and create to individual gene-level files
        :param indexes: (List) start/stop index points to read from/to
        :param matrix_cells: (List) cell names from matrix file
        :param cluster_cells: (List) cell names from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        start_pos, end_pos = indexes
        print(f"reading {self.matrix_file_path} at index {start_pos}:{end_pos}")
        with ExpressionWriter.open_file(self.matrix_file_path) as matrix_file:
            current_pos = start_pos
            matrix_file.seek(current_pos)
            while current_pos < end_pos:
                line = matrix_file.readline()
                self.process_dense_line(line, matrix_cells, cluster_cells, cluster_name, data_dir)
                current_pos += len(line)

    def process_dense_line(self, line, matrix_cells, cluster_cells, cluster_name, data_dir):
        """
        Process a single line from a dense matrix and write gene-level data

        :param matrix_cells: (List) cells from header line of dense matrix
        :param cluster_cells: (List) cell names from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir (String) name out output dir
        :param line: (String) single line from dense matrix
        """
        clean_line = line.rstrip().replace('"','')
        raw_vals = re.split(COMMA_OR_TAB, clean_line)
        gene_name = raw_vals.pop(0)
        exp_vals = [round(float(val), self.precision) if float(val) != 0.0 else 0 for val in raw_vals]
        filtered_expression = ExpressionWriter.filter_expression_for_cluster(
            cluster_cells, matrix_cells, exp_vals
        )
        ExpressionWriter.write_gene_scores(gene_name, cluster_name, filtered_expression, data_dir)

    @staticmethod
    def filter_expression_for_cluster(cluster_cells, exp_cells, exp_scores):
        """
        Assemble a List of expression scores, filtered & ordered from a List of cluster cells
        Will substitute 0 as a value for any cell not seen in the expression file

        :param cluster_cells: (List) cluster cell names
        :param exp_cells: (List) expression cell names
        :param exp_scores: (List) expression values, in the same order as exp_cells
        :return: (List) Expression values, ordered by cluster_cells
        """
        observed_exp = dict(zip(exp_cells, exp_scores))
        return (observed_exp.get(cell, 0) for cell in cluster_cells)

    @staticmethod
    def write_gene_scores(gene_name, cluster_name, exp_values, data_dir):
        """
        Write a JSON array of expression values
        Filename uses {cluster_name}--{gene_name}.json format

        :param gene_name: (String) Name of gene
        :param cluster_name: (String) Name of cluster
        :param exp_values: (List) expression values
        :param data_dir: (String) name out output dir
        """
        with gzip.open(f"{data_dir}/{cluster_name}--{gene_name}.json.gz", "wt") as file:
            json.dump(list(exp_values), file, separators=(',', ':'))


    def main(self):
        """
        Main handler, determines type of processing to execute (dense vs. sparse)
        """
        cluster_file_path = args.cluster_file
        expression_file_path = args.matrix_file
        sanitized_cluster_name = re.sub(r'\W', '_', args.cluster_name)
        data_dir = ExpressionWriter.make_data_dir(sanitized_cluster_name)
        cluster_cells = ExpressionWriter.get_cluster_cells(cluster_file_path)
        if self.genes_file is not None and self.barcodes_file is not None:
            print(f"reading {expression_file_path} as sparse matrix")
            genes_file = ExpressionWriter.open_file(args.genes_file)
            genes = ExpressionWriter.load_entities_as_list(genes_file, 1)
            barcodes_file = ExpressionWriter.open_file(args.barcodes_file)
            barcodes = ExpressionWriter.load_entities_as_list(barcodes_file)
            self.divide_sparse_matrix(genes, data_dir)
            self.write_empty_sparse_genes(genes, len(cluster_cells), sanitized_cluster_name, data_dir)
            self.process_sparse_data_fragments(barcodes, cluster_cells, sanitized_cluster_name, data_dir)
        else:
            print(f"reading {expression_file_path} as dense matrix")
            self.process_dense_data(cluster_cells, sanitized_cluster_name, data_dir)

if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster-file', help='path to cluster file', required=True)
    parser.add_argument('--cluster-name', help='name of cluster object', required=True)
    parser.add_argument('--matrix-file', help='path to matrix file', required=True)
    parser.add_argument('--precision', help='number of digits of precision for non-zero data',
                        nargs='?', const=1, type=int)
    parser.add_argument('--cores', help='number of cores to use in multiprocessing',
                        nargs='?', const=1, type=int)
    parser.add_argument('--genes-file', help='path to genes file (None for dense matrix files)')
    parser.add_argument('--barcodes-file', help='path to barcodes file (None for dense matrix files)')
    args = parser.parse_args()
    expression_file = args.matrix_file
    cluster_file = args.cluster_file
    writer = ExpressionWriter(
        args.matrix_file,
        args.cluster_file,
        args.cluster_name,
        args.genes_file,
        args.barcodes_file,
        args.precision,
        args.cores
    )
    writer.main()
    end_time = time.time()
    time_in_min = round(float(end_time - start_time), 3) / 60
    print(f"completed, total runtime in minutes: {time_in_min}")