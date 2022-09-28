"""expression_writer.py
extract gene-level expression data from dense/sparse matrix files and process in the
context of a cluster file's cells.
this mimics the `expression` data array in visualization API responses

EXAMPLES

dense matrix:
python3 expression_writer.py --matrix-file ../../tests/data/dense_expression_matrix.txt \
                             --cluster-file ../../tests/data/cluster_example.txt \
                             --cluster-name 'Dense Example'

sparse matrix:
python3 expression_writer.py --matrix-file ../../tests/data/mtx/matrix_with_header.mtx \
                             --genes-file ../../tests/data/mtx/sampled_genes.tsv \
                             --barcodes-file ../../tests/data/mtx/barcodes.tsv \
                             --cluster-file ../../tests/data/mtx/cluster_mtx_barcodes.tsv \
                             --cluster-name 'Sparse Example'
"""

import os
import re
import argparse
import time
import multiprocessing
import sys
from functools import partial

try:
    from writer_functions import round_exp, open_file, make_data_dir, get_matrix_size, get_cluster_cells,\
        load_entities_as_list, read_sparse_matrix_slice, process_sparse_fragment, write_gene_scores,\
        process_dense_line, COMMA_OR_TAB
    from monitor import setup_logger
except ImportError:
    from .writer_functions import round_exp, open_file, make_data_dir, get_matrix_size, get_cluster_cells, \
        load_entities_as_list, read_sparse_matrix_slice, process_sparse_fragment, write_gene_scores, \
        process_dense_line, COMMA_OR_TAB
    from .monitor import setup_logger

class ExpressionWriter:
    denominator = 2 if re.match('darwin', sys.platform) else 1
    num_cores = int(multiprocessing.cpu_count() / denominator) - 1

    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")

    def __init__(
        self, matrix_file_path, matrix_file_type, cluster_file_path, cluster_name, gene_file, barcode_file, **kwargs
    ):
        self.matrix_file_path = matrix_file_path
        self.matrix_file_type = matrix_file_type
        self.cluster_file_path = cluster_file_path
        self.cluster_name = cluster_name
        self.gene_file = gene_file
        self.barcode_file = barcode_file

    def get_file_seek_points(self):
        """
        Determine start/stop points in a matrix to process in parallel
        Will read in chunks and return a list of start/stop points
        Ensures breaks on newlines

        :return: (List<List>)
        """
        file_size = get_matrix_size(self.matrix_file_path)
        chunk_size = int(file_size / self.num_cores)
        ExpressionWriter.dev_logger.info(
            f" determining seek points for {self.matrix_file_path} with chunk size {chunk_size}"
        )
        with open_file(self.matrix_file_path) as matrix_file:
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
        ExpressionWriter.dev_logger.info(f" loading sparse data from {self.matrix_file_path}")
        slice_indexes = self.get_file_seek_points()
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(read_sparse_matrix_slice,
                            matrix_file_path=self.matrix_file_path, genes=genes, data_dir=data_dir)
        pool.map(processor, slice_indexes)


    def process_sparse_data_fragments(self, barcodes, cluster_cells, cluster_name, data_dir):
        """
        Find and process all generated single-gene sparse data fragments

        :param barcodes: (List) list of cell barcodes
        :param cluster_cells: (List) list of cells from cluster file
        :param cluster_name: (String) name of cluster object
        :param data_dir: (String) name out output dir
        """
        fragments = os.listdir(f"{data_dir}/gene_entries")
        ExpressionWriter.dev_logger.info(f" subdivision complete, processing {len(fragments)} fragments")
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(process_sparse_fragment,
                            barcodes=barcodes, cluster_cells=cluster_cells,
                            cluster_name=cluster_name, data_dir=data_dir)
        pool.map(processor, fragments)

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
        processor = partial(write_gene_scores,
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
        matrix_file = open_file(self.matrix_file_path)
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
        ExpressionWriter.dev_logger.info(f" reading {self.matrix_file_path} at index {start_pos}:{end_pos}")
        with open_file(self.matrix_file_path) as matrix_file:
            current_pos = start_pos
            matrix_file.seek(current_pos)
            while current_pos < end_pos:
                line = matrix_file.readline()
                process_dense_line(line, matrix_cells, cluster_cells, cluster_name, data_dir)
                current_pos += len(line)

    def render_artifacts(self):
        """
        Main handler, determines type of processing to execute (dense vs. sparse)
        """
        start_time = time.time()
        sanitized_cluster_name = re.sub(r'\W', '_', self.cluster_name)
        data_dir = make_data_dir(sanitized_cluster_name)
        cluster_cells = get_cluster_cells(self.cluster_file_path)
        if self.matrix_file_type == 'mtx' and self.gene_file is not None and self.barcode_file is not None:
            ExpressionWriter.dev_logger.info(f" reading {self.matrix_file_path} as sparse matrix")
            genes_file = open_file(args.gene_file)
            genes = load_entities_as_list(genes_file)
            barcodes_file = open_file(args.barcode_file)
            barcodes = load_entities_as_list(barcodes_file)
            self.divide_sparse_matrix(genes, data_dir)
            self.write_empty_sparse_genes(genes, len(cluster_cells), sanitized_cluster_name, data_dir)
            self.process_sparse_data_fragments(barcodes, cluster_cells, sanitized_cluster_name, data_dir)
        elif self.matrix_file_type == 'dense':
            ExpressionWriter.dev_logger.info(f" reading {self.matrix_file_path} as dense matrix")
            self.process_dense_data(cluster_cells, sanitized_cluster_name, data_dir)
        end_time = time.time()
        total_time = end_time - start_time
        time_in_min = round_exp(total_time, 3) / 60
        ExpressionWriter.dev_logger.info(f" completed, total runtime in minutes: {time_in_min}")

if __name__ == '__main__':
    """
    Main function to use as standalone module outside of ingest_pipeline.py
    """
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster-file', help='path to cluster file', required=True)
    parser.add_argument('--cluster-name', help='name of cluster object', required=True)
    parser.add_argument('--matrix-file-path', help='path to matrix file', required=True)
    parser.add_argument('--matrix-file-type', help='type to matrix file (dense or mtx)', required=True)
    parser.add_argument('--gene-file', help='path to gene file (None for dense matrix files)')
    parser.add_argument('--barcode-file', help='path to barcode file (None for dense matrix files)')
    args = parser.parse_args()
    expression_file = args.matrix_file
    cluster_file = args.cluster_file
    writer = ExpressionWriter(
        args.matrix_file, args.matrix_file_type, args.cluster_file, args.cluster_name, args.genes_file,
        args.barcodes_file
    )
    writer.render_artifacts()
    end_time = time.time()
    total_time = end_time - start_time
    time_in_min = round_exp(total_time, 3) / 60
    ExpressionWriter.dev_logger.info(f" completed, total runtime in minutes: {time_in_min}")