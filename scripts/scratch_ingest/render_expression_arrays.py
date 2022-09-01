#! /usr/bin/env python3

"""render_expression_arrays.py
extract gene-level expression data from dense/sparse matrix files and process in the
context of a cluster file's cells.
this mimics the `expression` data array in visualization API responses

EXAMPLES

dense matrix:
python3 render_expression_arrays.py --matrix-file ../../tests/data/dense_expression_matrix.txt \
                                    --cluster-file ../../tests/data/cluster_example.txt \
                                    --cluster-name 'Dense Example'

dense matrix with precision override:
python3 render_expression_arrays.py --matrix-file ../../tests/data/dense_expression_matrix.txt \
                                    --cluster-file ../../tests/data/cluster_example.txt \
                                    --cluster-name 'Dense Example' --precision 1

sparse matrix:
python3 render_expression_arrays.py --matrix-file ../../tests/data/mtx/sorted_matrix_header.mtx \
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
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# regex to split on commas and tabs
COMMA_OR_TAB = r"[,\t]"

# regex to split on comma/tab/space
ALL_DELIM = r"[\s,\t]"

# Gzip magic number
GZIP_MAGIC_NUM = b'\x1f\x8b'

# default level of precision
precision = 3

num_cores = multiprocessing.cpu_count() - 1

def is_gz_file(filepath):
    """
    Determine if a file is gzipped by reading the first two bytes and comparing to the 'magic number'
    Args:
        filepath (String): path to file to check

    Returns:
        (Bool): T/F if file is gzipped
    """
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == GZIP_MAGIC_NUM
        
def open_file(filepath):
    """
    Open a file with the correct reader
    Args:
        filepath (String): path to file to open

    Returns:
        (TextIOWrapper): opened file
    """
    if is_gz_file(filepath):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def make_data_dir(name):
    """
    Make a directory to put output files in.  Appends a UUID to the end of the directory to allow for
    quick iteration and side-by-side comparison of outputs without manually moving directories
    Args:
        name (String): name of directory
        
    Returns:
        (String): name of created diretory with UUID appended to the end
    """
    dirname = str(f"{name}-{uuid.uuid4()}")
    print(f"creating data directory at {dirname}")
    os.mkdir(dirname)
    os.mkdir(f"{dirname}/gene_entries")
    return dirname
    
def get_cluster_cells(file_path):
    """
    Return a list of cell names from a cluster file (1st column, starting line 3)
    Args:
        file_path (String): absolute path to cluster_file
        
    Returns:
        (List): cluster cell names
    """
    cells = []
    with open_file(file_path) as cluster_file:
        cluster_file.readline()
        cluster_file.readline()
        for row in cluster_file:
            cell = re.split(COMMA_OR_TAB, row)[0]
            cells.append(cell)
    return cells

def load_entities_as_list(file, column = None):
    """
    Read an entire 10X feature/barcode file into a list for parsing sparse data
    Args:
        file (TextIOWrapper): open file object
        
    Returns:
        (List): entities from file
    """
    print(f"reading entities from {file.name}")
    if not column:
        return list(map(str.rstrip, file.readlines()))
    else:
        return list(re.split(ALL_DELIM, line.strip())[column] for line in file)

def divide_sparse_matrix(matrix_file_path, genes, data_dir):
    print(f"loading sparse data from {matrix_file_path}")
    matrix_file = open_file(matrix_file_path)
    matrix_file.readline()
    matrix_file.readline()
    matrix_file.readline()
    pool = multiprocessing.Pool(num_cores)
    file_iterator = enumerate(matrix_file)
    processor = partial(process_sparse_line, genes, data_dir)
    pool.map(processor, file_iterator)

def process_sparse_line(genes, data_dir, processor_tuple):
    line_no = processor_tuple[0]
    line = processor_tuple[1]
    if line_no % 1000 == 0:
        print(f"at line {line_no}")
    gene_idx = int(line.split()[0])
    gene_name = genes[gene_idx - 1]
    outfile = f"{data_dir}/gene_entries/{gene_name}__entries.txt"
    write_data_to_fragment(line, outfile)

def write_data_to_fragment(fragment, file_path):
    with open(file_path, 'a+') as file:
        file.write(f"{fragment}")

def process_fragment(barcodes, cluster_cells, cluster_name, data_dir, fragment_name):
    """
    Process a single-gene sparse matrix fragment and write expression array
    Args:
        barcodes (List): list of cell barcodes
        cluster_cells (List): list of cells from cluster file
        cluster_name (String): name of cluster object
        data_dir (String): output path to write data in
        fragment_name (String): name of gene-level fragment file
    """
    gene_name = fragment_name.split('__')[0]
    full_path = f"{data_dir}/gene_entries/{fragment_name}"
    observed_cells = []
    exp_vals = []
    # write array of 0 values if file is empty
    if os.stat(full_path).st_size == 0:
        empty_exp = [0] * len(cluster_cells)
        write_gene_scores(cluster_name, gene_name, empty_exp, data_dir)
    for line in open_file(full_path):
        gene_idx, barcode_idx, exp_val = extract_sparse_line(line)
        observed_cells.append(barcodes[barcode_idx - 1])
        exp_vals.append(exp_val)
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, observed_cells, exp_vals
    )
    write_gene_scores(cluster_name, gene_name, filtered_expression, data_dir)

def process_sparse_data_fragments(barcodes, cluster_cells, cluster_name, data_dir):
    """
    Find and process all generated single-gene sparse data fragments
    Args:
        barcodes (List): list of cell barcodes
        cluster_cells (List): list of cells from cluster file
        cluster_name (String): name of cluster object
        data_dir (String): output path to write data in
    """
    fragments = os.listdir(f"{data_dir}/gene_entries")
    print(f"subdivision complete, processing {len(fragments)} fragments")
    pool = multiprocessing.Pool(num_cores)
    processor = partial(process_fragment, barcodes, cluster_cells, cluster_name, data_dir)
    pool.map(processor, fragments)

def extract_sparse_line(line):
    """
    Process a single line from a sparse matrix and extract values as integers
    Args:
        line (String): single line from matrix file
        
    Returns:
        (List): values as integers
    """
    gene_idx, barcode_idx, raw_exp = line.rstrip().split(' ')
    return [int(gene_idx), int(barcode_idx), round(float(raw_exp), precision)]

def process_dense_data(matrix_file_path, cluster_cells, cluster_name, data_dir):
    """
    Main handler to read dense matrix data and process entries at the gene level
    Args:
        matrix_file_path (String): path to dense matrix file
        cluster_cells (List): cell names from cluster file
        cluster_name (String): name of cluster object
        data_dir (String): output data directory
    """
    pool = multiprocessing.Pool(num_cores)
    with open_file(matrix_file_path) as matrix_file:
        header = matrix_file.readline().rstrip()
        values = re.split(COMMA_OR_TAB, header)
        cells = values[1:]
        processor = partial(process_dense_line, cells, cluster_cells, cluster_name, data_dir)
        pool.map(processor, (matrix_file.readlines))

def process_dense_line(matrix_cells, cluster_cells, cluster_name, data_dir, line):
    clean_line = line.rstrip()
    raw_vals = re.split(COMMA_OR_TAB, clean_line)
    gene_name = raw_vals.pop(0)
    exp_vals = [round(float(val), precision) if float(val) != 0.0 else 0 for val in raw_vals]
    filtered_expression = filter_expression_for_cluster(
        cluster_cells, matrix_cells, exp_vals
    )
    write_gene_scores(cluster_name, gene_name, filtered_expression, data_dir)

def filter_expression_for_cluster(cluster_cells, exp_cells, exp_scores):
    """
    Assemble a List of expression scores, filtered & ordered from a List of cluster cells
    Will substitute 0 as a value for any cell not seen in the expression file
    
    Args:
        cluster_cells (List): cluster cell names
        exp_cells (List): expression cell names
        exp_scores (List): expression values, in the same order as exp_cells
        
    Returns:
        (List): Expression values, ordered by cluster_cells
    """
    observed_exp = dict(zip(exp_cells, exp_scores))
    return (observed_exp.get(cell, 0) for cell in cluster_cells)


def write_gene_scores(cluster_name, gene_name, exp_values, data_dir):
    """
    Write a JSON array of expression values
    Filename uses {cluster_name}--{gene_name}.json format

    Args:
        cluster_name (String): Name of cluster
        gene_name (String): Name of gene
        exp_values (List): expression values
        data_dir (String): output data directory to write in
        
    Returns:
        (File): JSON file as a List of expression values
    """
    with gzip.open(f"{data_dir}/{cluster_name}--{gene_name}.json.gz", "wt") as file:
        json.dump(list(exp_values), file, separators=(',', ':'))

if __name__ == '__main__':
    # parse cli arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster-file', help='path to cluster file', required=True)
    parser.add_argument('--cluster-name', help='name of cluster object', required=True)
    parser.add_argument('--matrix-file', help='path to matrix file', required=True)
    parser.add_argument('--precision', help='number of digits of precision for non-zero data')
    parser.add_argument('--genes-file', help='path to genes file (None for dense matrix files)')
    parser.add_argument('--barcodes-file', help='path to barcodes file (None for dense matrix files)')
    args = parser.parse_args()

    # main execution, set cluster name
    start_time = time.time()
    cluster_file_path = args.cluster_file
    expression_file_path = args.matrix_file
    sanitized_cluster_name = re.sub(r'\W', '_', args.cluster_name)
    data_dir = make_data_dir(sanitized_cluster_name)
    cluster_cells = get_cluster_cells(cluster_file_path)
    print(f"using {num_cores} cores for multiprocessing")
    if args.precision is not None:
        precision = int(args.precision)
        print(f"using {precision} digits of precision for non-zero data")
    if args.genes_file is not None and args.barcodes_file is not None:
        print(f"reading {expression_file_path} as sparse matrix")
        genes_file = open_file(args.genes_file)
        genes = load_entities_as_list(genes_file, 1)
        barcodes_file = open_file(args.barcodes_file)
        barcodes = load_entities_as_list(barcodes_file)
        divide_sparse_matrix(expression_file_path, genes, data_dir)
        process_sparse_data_fragments(barcodes, cluster_cells, sanitized_cluster_name, data_dir)
    else:
        print(f"reading {expression_file_path} as dense matrix")
        process_dense_data(
            expression_file_path, cluster_cells, sanitized_cluster_name, data_dir
        )

    end_time = time.time()
    time_in_min = round(float(end_time - start_time), 3) / 60
    print(f"completed, total runtime in minutes: {time_in_min}")
