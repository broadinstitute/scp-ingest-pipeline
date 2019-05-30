"""Command-line interface for ingesting MTX file bundles into Firestore

DESCRIPTION
This CLI extracts and transforms gene expression data in an MTX file bundle,
and loads it to Google Cloud Firestore.  An MTX file bundle consists of A) an 
.mtx file in Matrix Market matrix coordinate format, B) a genes.tsv file,
and a barcodes.tsv file.  These commonly provided from 10x Genomics v2.

PREREQUISITES
You must have google cloud firestore installed, authenticated
and configured.

EXAMPLES
# Transform and load MTX matrix bundle
$ python mtx_ingest.py --input-dir ../tests/data --matrix-bundle matrix.mtx genes.tsv barcodes.tsv
"""

import argparse
import sys
import json
import os
import time
from itertools import islice

from google.cloud import firestore
import scipy.io

db = firestore.Client()

parser = argparse.ArgumentParser(
    prog = 'mtx_ingest.py'
)

parser.add_argument(
    '--input-dir', default='.',
    help='Input directory for matrix bundle files'
)
parser.add_argument(
    '--matrix-bundle', default=None, nargs='+',
    help='Names of .mtx, .genes.tsv, and .barcodes.tsv files'
)

args = parser.parse_args()
input_dir = args.input_dir
matrix_bundle = args.matrix_bundle

def parse_mtx_bundle(bundle_paths):
    """Get relevant iterables for each file of the MTX bundle
    """
    mtx_path, genes_path, barcodes_path = bundle_paths
    with open(mtx_path) as f:
        matrix = scipy.io.mmread(mtx_path)
    with open(genes_path) as f:
        genes = [g.strip() for g in f.readlines()]
    with open(barcodes_path) as f:
        cells = [c.strip() for c in f.readlines()]
    return matrix, genes, cells

def map_genes_to_expression_values(bundle_paths):    
    matrix, genes, cells = parse_mtx_bundle(bundle_paths)
    exp_by_gene = {}
    for i, j, exp_score in zip(matrix.row, matrix.col, matrix.data):
        gene = genes[i]
        exp_score = float(exp_score)
        if gene in exp_by_gene:
            exp_by_gene[gene].append(exp_score)
        else:
            exp_by_gene[gene] = [exp_score]
    return cells, exp_by_gene

    
## Save to firestore
## Abstract this out
def add_data_to_firestore(exp_by_gene, cells, matrix_path):
    batch = db.batch()
    for gene, exp_vals in exp_by_gene.items():
        doc_ref = db.collection('genes').document(gene)
        doc_data = {
            'source_file_name': matrix_path,
            'data_arrays': {
                'cell_names': cells,
                'expression_scores': exp_vals
            }
        }
        batch.set(doc_ref, doc_data)
    batch.commit()
    time.sleep(2)

bundle_paths = [os.path.join(input_dir, fn) for fn in matrix_bundle]

cells, exp_by_gene = map_genes_to_expression_values(bundle_paths)
add_data_to_firestore(exp_by_gene, cells, matrix_bundle)