"""expression_writer.py
extract gene-level expression data from dense/sparse matrix files and process in the
context of a cluster file's cells.
this mimics the `expression` data array in visualization API responses

EXAMPLES (must be invoked via ingest_pipeline.py)

dense matrix:
python3 ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 \
                             render_expression_arrays --matrix-file-path ../tests/data/dense_expression_matrix.txt \
                             --matrix-file-type dense \
                             --cluster-file ../tests/data/cluster_example.txt \
                             --cluster-name 'Dense Example' --render-expression-arrays

sparse matrix:
python3 ingest_pipeline.py --study-id 5d276a50421aa9117c982845 --study-file-id 5dd5ae25421aa910a723a337 \
                             render_expression_arrays --matrix-file-path ../tests/data/mtx/matrix_with_header.mtx \
                             --matrix-file-type mtx \
                             --gene-file ../tests/data/mtx/sampled_genes.tsv \
                             --barcode-file ../tests/data/mtx/barcodes.tsv \
                             --cluster-file ../tests/data/mtx/cluster_mtx_barcodes.tsv \
                             --cluster-name 'Sparse Example' --render-expression-arrays
"""

import os
import subprocess
import re
import multiprocessing
import sys
import datetime
from dateutil.relativedelta import relativedelta
from functools import partial

try:
    from writer_functions import round_exp, encode_cluster_name, open_file, make_data_dir, get_matrix_size, \
        get_cluster_cells, load_entities_as_list, process_sparse_fragment, write_gene_scores, process_dense_line, \
        COMMA_OR_TAB
    from monitor import setup_logger
    from ingest_files import IngestFiles
except ImportError:
    from .writer_functions import round_exp, encode_cluster_name, open_file, make_data_dir, get_matrix_size, \
        get_cluster_cells, load_entities_as_list, process_sparse_fragment, write_gene_scores, process_dense_line, \
        COMMA_OR_TAB
    from .monitor import setup_logger
    from .ingest_files import IngestFiles

class ExpressionWriter:
    DELOCALIZE_FOLDER = "_scp_internal/cache/expression_scatter/data"
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

        # localize main files, if needed
        self.local_matrix_path = open_file(self.matrix_file_path)[1]
        self.local_cluster_path = open_file(self.cluster_file_path)[1]

        # set storage bucket name, if needed
        self.bucket = self.get_storage_bucket_name()

    def get_storage_bucket_name(self):
        """
        Load GCS storage bucket, if available
        :return: (google.cloud.storage.Bucket)
        """
        path_header = self.matrix_file_path[:5]
        if path_header == 'gs://':
            path_segments = self.matrix_file_path[5:].split("/")
            bucket_name = path_segments[0]
            return f"{path_header}{bucket_name}"

    def get_file_seek_points(self):
        """
        Determine start/stop points in a matrix to process in parallel
        Will read in chunks and return a list of start/stop points
        Ensures breaks on newlines

        :return: (List<List>)
        """
        file_size = get_matrix_size(self.local_matrix_path)
        chunk_size = int(file_size / self.num_cores)
        self.dev_logger.info(
            f" determining seek points for {self.local_matrix_path} with chunk size {chunk_size}"
        )
        with open_file(self.local_matrix_path)[0] as matrix_file:
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
        :param data_dir: (String) name of output dir
        """
        self.dev_logger.info(f" loading sparse data from {self.local_matrix_path}")
        slice_indexes = self.get_file_seek_points()
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(self.read_sparse_matrix_slice, genes=genes, data_dir=data_dir)
        pool.map(processor, slice_indexes)

    def read_sparse_matrix_slice(self, indexes, genes, data_dir):
        """
        Read a sparse matrix using start/stop indexes and append to individual gene-level files

        :param indexes: (List) start/stop index points to read from/to
        :param genes: (List) gene names from features file
        :param data_dir: (String) name of output dir
        """
        start_pos, end_pos = indexes
        self.dev_logger.info(f"reading {self.local_matrix_path} at index {start_pos}:{end_pos}")
        with open_file(self.local_matrix_path)[0] as matrix_file:
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

    def process_sparse_data_fragments(self, barcodes, cluster_cells, data_dir):
        """
        Find and process all generated single-gene sparse data fragments

        :param barcodes: (List) list of cell barcodes
        :param cluster_cells: (List) list of cells from cluster file
        :param data_dir: (String) name of output dir
        """
        fragments = os.listdir(f"{data_dir}/gene_entries")
        self.dev_logger.info(f" subdivision complete, processing {len(fragments)} fragments")
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(process_sparse_fragment, barcodes=barcodes, cluster_cells=cluster_cells, data_dir=data_dir)
        pool.map(processor, fragments)

    def write_empty_sparse_genes(self, genes, num_cluster_cells, data_dir):
        """
        Write out empty arrays of expression values for genes with no significant expression in a sparse matrix

        :param genes: (List) gene names from features file
        :param num_cluster_cells: (Integer) number of cells from cluster file
        :param data_dir: (String) name of output dir
        """
        gene_fragments = filter(lambda file: file[0] != '.', os.listdir(f"{data_dir}/gene_entries"))
        significant_genes = set([gene.split('__')[0] for gene in gene_fragments])
        missing_genes = [gene for gene in genes if gene not in significant_genes]
        empty_expression = [0] * num_cluster_cells
        pool = multiprocessing.Pool(self.num_cores)
        processor = partial(write_gene_scores, exp_values=empty_expression, data_dir=data_dir)
        pool.map(processor, missing_genes)

    def process_dense_data(self, cluster_cells, data_dir):
        """
        Main handler to read dense matrix data and process entries at the gene level

        :param cluster_cells: (List) cell names from cluster file
        :param data_dir: (String) name of output dir
        """
        pool = multiprocessing.Pool(self.num_cores)
        slice_indexes = self.get_file_seek_points()
        matrix_file, local_matrix_path = open_file(self.matrix_file_path)
        header = matrix_file.readline().rstrip()
        values = re.split(COMMA_OR_TAB, header)
        cells = values[1:]
        processor = partial(
            self.read_dense_matrix_slice, matrix_cells=cells, cluster_cells=cluster_cells, data_dir=data_dir
        )
        pool.map(processor, slice_indexes)

    def read_dense_matrix_slice(self, indexes, matrix_cells, cluster_cells, data_dir):
        """
        Read a dense matrix using start/stop indexes and create to individual gene-level files
        :param indexes: (List) start/stop index points to read from/to
        :param matrix_cells: (List) cell names from matrix file
        :param cluster_cells: (List) cell names from cluster file
        :param data_dir: (String) name of output dir
        """
        start_pos, end_pos = indexes
        self.dev_logger.info(f" reading {self.local_matrix_path} at index {start_pos}:{end_pos}")
        with open_file(self.local_matrix_path)[0] as matrix_file:
            current_pos = start_pos
            matrix_file.seek(current_pos)
            while current_pos < end_pos:
                line = matrix_file.readline()
                process_dense_line(line, matrix_cells, cluster_cells, data_dir)
                current_pos += len(line)

    def render_artifacts(self):
        """
        Main handler, determines type of processing to execute (dense vs. sparse)
        """
        start_time = datetime.datetime.now()
        sanitized_cluster_name = encode_cluster_name(self.cluster_name)
        make_data_dir(sanitized_cluster_name)
        cluster_cells = get_cluster_cells(self.local_cluster_path)
        if self.matrix_file_type == 'mtx' and self.gene_file is not None and self.barcode_file is not None:
            self.dev_logger.info(f" reading {self.matrix_file_path} as sparse matrix")
            genes_file = open_file(self.gene_file)[0]
            genes = load_entities_as_list(genes_file)
            barcodes_file = open_file(self.barcode_file)[0]
            barcodes = load_entities_as_list(barcodes_file)
            self.divide_sparse_matrix(genes, sanitized_cluster_name)
            self.write_empty_sparse_genes(genes, len(cluster_cells), sanitized_cluster_name)
            self.process_sparse_data_fragments(barcodes, cluster_cells, sanitized_cluster_name)
        elif self.matrix_file_type == 'dense':
            self.dev_logger.info(f" reading {self.matrix_file_path} as dense matrix")
            self.process_dense_data(cluster_cells, sanitized_cluster_name)
        end_time = datetime.datetime.now()
        time_diff = relativedelta(end_time, start_time)
        self.dev_logger.info(
            f" completed, total runtime: {time_diff.hours}h, {time_diff.minutes}m, {time_diff.seconds}s"
        )
        self.delocalize_outputs(sanitized_cluster_name)

    def delocalize_outputs(self, cluster_name):
        """
        Copy all output files to study bucket in parallel using gsutil (since there are usually ~25-30K files)

        :param data_dir: (String) name of output dir
        :param cluster_name: (String) encoded name of cluster
        """
        if self.bucket is not None:
            bucket_path = f"{self.bucket}/{self.DELOCALIZE_FOLDER}/{cluster_name}"
            self.dev_logger.info(f" pushing all output files to {bucket_path}")
            subprocess.Popen(
                ["gsutil", "-h" , "Content-Encoding:gzip", "-m", "cp", f"{cluster_name}/*.json.gz", f"{bucket_path}"]
            )
            self.dev_logger.info(" push completed")

