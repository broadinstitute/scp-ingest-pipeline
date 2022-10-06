from __future__ import annotations
from typing import TextIO
import re
import os
import gzip
import json

try:
    from ingest_files import IngestFiles
except ImportError:
    from .ingest_files import IngestFiles

# regex to split on commas and tabs
COMMA_OR_TAB = r"[,\t]"

# regex to split on comma/tab/space
ALL_DELIM = r"[\s,\t]"

# Gzip magic number
GZIP_MAGIC_NUM = b'\x1f\x8b'

ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

def is_gz_file(filepath) -> bool:
    """
    Determine if a file is gzipped by checking first two bytes against GZIP_MAGIC_NUM

    :param filepath: (str) path to file to check
    :returns: (bool)
    """
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == GZIP_MAGIC_NUM

def round_exp(val, precision) -> float:
    """
    Round a raw value to a specified precision

    :param val: (str) raw value to convert and round
    :param precision: (int) number of digits of precision
    :returns: (float)
    """
    return round(float(val), precision)

def encode_cluster_name(name) -> str:
    """
    Encodes a cluster name to be used in a HTTP request
    Will encode + signs as 'pos' to differentiate + and -
    :param name: (str) name of cluster
    :returns: (str) encoded cluster name
    """
    plus_converted_string = name.replace('+', 'pos')
    return re.sub(r'\W', '_', plus_converted_string)

def open_file(filepath) -> list[TextIO, str]:
    """
    Open a file with the correct reader, e.g. even if it's gzipped

    :param filepath: (str) path to file to open
    :returns: list[TextIO, str]
    """
    ingest_file = IngestFiles(filepath, ALLOWED_FILE_TYPES)
    file_io, local_path = ingest_file.resolve_path(filepath)
    return file_io, local_path

def make_data_dir(name):
    """
    Make a directory to put output/work files in

    :param name: (str) name of directory
    """
    os.mkdir(name)
    os.mkdir(f"{name}/gene_entries") # for slicing sparse matrix files

def get_matrix_size(matrix_file_path) -> int:
    """
    Get the actual size of a matrix file, whether gzipped or not
    Because gzipped files only use 4 bytes to store the 'original' file size,
    any files larger that 4 GiB will only store the modulus of 2 ** 32, not the
    actual size of the file.
    Using seek() and tell() will return the last byte position

    :param matrix_file_path: (str) path to matrix file
    :returns: (int)
    """
    matrix_file, local_path = open_file(matrix_file_path)
    if is_gz_file(local_path):
        matrix_file.seek(0, 2)
        return matrix_file.tell()
    else:
        return os.stat(local_path).st_size

def get_cluster_cells(file_path) -> list:
    """
    Return a list of cell names from a cluster file (1st column, starting line 3)

    :param file_path: (str) absolute path to cluster_file
    :returns: (list)
    """
    cells = []
    with open_file(file_path)[0] as file:
        file.readline()
        file.readline()
        for row in file:
            cell = re.split(COMMA_OR_TAB, row)[0]
            cells.append(cell)
    return cells

def load_entities_as_list(file) -> list:
    """
    Read an entire 10X feature/barcode file into a list for parsing sparse data

    :param file: (TextIO) open file object
    :param column: (int) specific column to extract from entity file
    """
    column = get_entity_index(file)
    return list(re.split(ALL_DELIM, line.strip())[column] for line in file)

def get_entity_index(file) -> int:
    """
    Determine which column from a 10X entity file contains valid data

    :param file: (TextIO) open file object
    """
    first_line = file.readline().strip()
    file.seek(0) # gotcha to reset pointer to beginning
    entities = re.split(ALL_DELIM, first_line)
    return 0 if len(entities) == 1 else 1

def process_sparse_fragment(fragment_name, barcodes, cluster_cells, data_dir):
    """
    Extract and filter a single-gene sparse matrix fragment and write expression array

    :param fragment_name: (str) name of gene-level fragment file
    :param barcodes: (list) list of cell barcodes
    :param cluster_cells: (list) list of cells from cluster file
    :param data_dir: (str) name out output dir
    """
    gene_name = fragment_name.split('__')[0]
    full_path = f"{data_dir}/gene_entries/{fragment_name}"
    observed_cells = []
    exp_vals = []
    for line in open(full_path):
        gene_idx, barcode_idx, exp_val = extract_sparse_line(line)
        observed_cells.append(barcodes[barcode_idx - 1])
        exp_vals.append(exp_val)
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, observed_cells, exp_vals
    )
    write_gene_scores(gene_name, filtered_expression, data_dir)

def extract_sparse_line(line) -> list:
    """
    Process a single line from a sparse matrix and extract values as integers

    :param line: (str) single line from matrix file
    :returns: values for gene index, barcode index, and rounded expression
    """
    gene_idx, barcode_idx, raw_exp = line.rstrip().split(' ')
    return [int(gene_idx), int(barcode_idx), round_exp(raw_exp, 3)]

def process_dense_line(line, matrix_cells, cluster_cells, data_dir):
    """
    Process a single line from a dense matrix and write gene-level data

    :param matrix_cells: (list) cells from header line of dense matrix
    :param cluster_cells: (list) cell names from cluster file
    :param data_dir (str) name out output dir
    :param line: (str) single line from dense matrix
    """
    clean_line = line.rstrip().replace('"', '')
    raw_vals = re.split(COMMA_OR_TAB, clean_line)
    gene_name = raw_vals.pop(0)
    exp_vals = [round_exp(val, 3) if float(val) != 0.0 else 0 for val in raw_vals]
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, matrix_cells, exp_vals
    )
    if gene_name:
        write_gene_scores(gene_name, filtered_expression, data_dir)

def filter_expression_for_cluster(cluster_cells, exp_cells, exp_scores) -> list:
    """
    Assemble a List of expression scores, filtered & ordered from a List of cluster cells
    Will substitute 0 as a value for any cell not seen in the expression file

    :param cluster_cells: (list) cluster cell names
    :param exp_cells: (list) expression cell names
    :param exp_scores: (list) expression values, in the same order as exp_cells
    :returns: Expression values, ordered by cluster_cells
    """
    observed_exp = dict(zip(exp_cells, exp_scores))
    return (observed_exp.get(cell, 0) for cell in cluster_cells)

def write_gene_scores(gene_name, exp_values, data_dir):
    """
    Write a JSON array of expression values

    :param gene_name: (str) Name of gene
    :param exp_values: (list) expression values
    :param data_dir: (str) name out output dir
    """
    with gzip.open(f"{data_dir}/{gene_name}.json", "wt") as file:
        json.dump(list(exp_values), file, separators=(',', ':'))
