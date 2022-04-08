#! /usr/bin/env python
"""Sort sparse expression matrix files by gene, barcode
"""

import pandas as pd
import os
import sys
import argparse
import gzip
import shutil

def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, block_size)


def sort_sparse_matrix(matrix_file, sorted_matrix_file=None):
    """Sort the sparse matrix file by gene, barcode
        Arguments:
        matrix_file- sparse matrix file to sort
        sorted_matrix_file- optional output path of sorted file
        Outputs-
        output_name- gene, barcode sorted sparse matrix file, name is either sorted_matrix_file if provided, or "gene_sorted-" + matrix_file
    """
    if ('.mtx.gz' in matrix_file):
        unzipped_file = matrix_file.split('.gz')[0]
        print('Unzipping matrix file')
        gunzip_shutil(matrix_file, unzipped_file)
        matrix_file = unzipped_file

    if sorted_matrix_file:
        # if the sorted path name is provided, use it as the output name
        output_name = sorted_matrix_file
    else:
        # otherwise if the sorted path name is not provided, add "gene_sorted-" to the original path name
        directory, base_name = os.path.split(matrix_file)
        output_name = os.path.join(directory, "gene_sorted-" + base_name)
    # read sparse matrix
    print("Reading Sparse Matrix")
    headers = []
    with open(matrix_file) as matrix:
        line = next(matrix)
        while line.startswith("%"):
            headers = headers + [line]
            line = next(matrix)
        headers = headers + [line]
        df = pd.read_table(matrix, sep="\s+", names=['genes', 'barcodes', 'expr'])
    # sort sparse matrix
    print("Sorting Sparse Matrix")
    df = df.sort_values(by=['genes', 'barcodes'])
    # save sparse matrix
    print("Saving Sparse Matrix to:", output_name)
    with open(output_name, "w+") as output:
        output.write(''.join(headers))
    df.to_csv(output_name, sep=' ', index=False, header=0, mode="a")


def __main__(argv):
    """Command Line parser for sort_sparse_matrix
    Inputs-
    command line arguments
    """
    # create the argument parser
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # add arguments
    parser.add_argument('matrix_file', help='Sparse Matrix file')
    parser.add_argument('--sorted_matrix_file', '-o', help='Gene sorted sparse matrix file path', default=None)
    # call sort_sparse_matrix with parsed args
    args = parser.parse_args()
    sort_sparse_matrix(matrix_file=args.matrix_file, sorted_matrix_file=args.sorted_matrix_file)

# python default
if __name__ == '__main__':
    __main__(sys.argv)