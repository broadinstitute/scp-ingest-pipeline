import re
import os
import uuid
import gzip
import json

# regex to split on commas and tabs
COMMA_OR_TAB = r"[,\t]"

# regex to split on comma/tab/space
ALL_DELIM = r"[\s,\t]"

# Gzip magic number
GZIP_MAGIC_NUM = b'\x1f\x8b'

def is_gz_file(filepath):
    """
    Determine if a file is gzipped by checking first two bytes against GZIP_MAGIC_NUM

    :param filepath: (String) path to file to check
    :returns (Bool)
    """
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == GZIP_MAGIC_NUM

def round_exp(val, precision):
    """
    Round a raw value to a specified precision

    :param val: (String) raw value to convert and round
    :param precision: (Integer) number of digits of precision
    :return: (Float)
    """
    return round(float(val), precision)

def open_file(filepath):
    """
    Open a file with the correct reader, e.g. even if it's gzipped

    :param filepath: (String) path to file to open
    :returns: (TextIOWrapper)
    """
    if is_gz_file(filepath):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

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

def get_matrix_size(matrix_file_path):
    """
    Get the actual size of a matrix file, whether gzipped or not
    Because gzipped files only use 4 bytes to store the 'original' file size,
    any files larger that 4 GiB will only store the modulus of 2 ** 32, not the
    actual size of the file.
    Using seek() and tell() will return the last byte position

    :param matrix_file_path: (String) path to matrix file
    :return: (Int)
    """

    if is_gz_file(matrix_file_path):
        with open_file(matrix_file_path) as matrix_file:
            matrix_file.seek(0, 2)
            return matrix_file.tell()
    else:
        return os.stat(matrix_file_path).st_size

def get_cluster_cells(file_path):
    """
    Return a list of cell names from a cluster file (1st column, starting line 3)

    :param file_path: (String) absolute path to cluster_file
    :returns: (List)
    """
    cells = []
    with open_file(file_path) as file:
        file.readline()
        file.readline()
        for row in file:
            cell = re.split(COMMA_OR_TAB, row)[0]
            cells.append(cell)
    return cells

def load_entities_as_list(file):
    """
    Read an entire 10X feature/barcode file into a list for parsing sparse data

    :param file: (TextIOWrapper) open file object
    :param column: (Integer) specific column to extract from entity file
    :returns: (List)
    """
    column = get_entity_index(file)
    print(f"reading entities from {file.name} in column {column + 1}")
    return list(re.split(ALL_DELIM, line.strip())[column] for line in file)

def get_entity_index(file):
    """
    Determine which column from a 10X entity file contains valid data

    :param file: (TextIOWrapper) open file object
    :return: (Integer)
    """
    first_line = file.readline().strip()
    file.seek(0) # gotcha to reset pointer to beginning
    entities = re.split(ALL_DELIM, first_line)
    return 0 if len(entities) == 1 else 1

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
    with open_file(matrix_file_path) as matrix_file:
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

def process_sparse_fragment(fragment_name, barcodes, cluster_cells, cluster_name, data_dir):
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
    for line in open_file(full_path):
        gene_idx, barcode_idx, exp_val = extract_sparse_line(line)
        observed_cells.append(barcodes[barcode_idx - 1])
        exp_vals.append(exp_val)
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, observed_cells, exp_vals
    )
    write_gene_scores(gene_name, cluster_name, filtered_expression, data_dir)

def extract_sparse_line(line):
    """
    Process a single line from a sparse matrix and extract values as integers

    :param line: (String) single line from matrix file
    :return: (List): values as integers
    """
    gene_idx, barcode_idx, raw_exp = line.rstrip().split(' ')
    return [int(gene_idx), int(barcode_idx), round_exp(raw_exp, 3)]

def process_dense_line(line, matrix_cells, cluster_cells, cluster_name, data_dir):
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
    exp_vals = [round_exp(val, 3) if float(val) != 0.0 else 0 for val in raw_vals]
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, matrix_cells, exp_vals
    )
    write_gene_scores(gene_name, cluster_name, filtered_expression, data_dir)

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
