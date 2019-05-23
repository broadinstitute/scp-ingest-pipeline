"""Command-line interface for ingesting dense matrix files into Firestore

DESCRIPTION
This CLI maps cell ids to expression values, and gene ids from
a zarr file and puts them into firestore.

PREREQUISITES
You must have google cloud firestore installed, authenticated
and configured.

EXAMPLES
# Transform dense matrix and load to Firestore
$ python dense_ingest.py  200k_subsample_4k_PMBC.zarr
"""

import argparse
import sys
import json
import os
import time
import logging
from itertools import islice

from google.cloud import firestore

start_time = time.time()

db = firestore.Client()

parser = argparse.ArgumentParser(
    prog = 'dense_ingest.py'
)

parser.add_argument(
    '--matrix-path', default=None,
    help='Path to dense matrix file'
)

import logging

def get_logger(output_dir, log_name):
    """Creates a log file and returns an object to interface with it.
    """
    logger = logging.getLogger(log_name)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(output_dir + log_name + '.log')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

args = parser.parse_args()

def map_genes_to_expression_values(file_path):    
    exp_by_gene = {} # Expression values by gene name
    with open(file_path) as f:
        lines = f.readlines()
    cell_names = lines[0].split('\t')[1:]
    for line in lines[1:]: # Skip header
        columns = line.strip().split('\t') # TODO: Abstract delimiter
        gene = columns[0]
        expression_values = columns[1:]
        exp_by_gene[gene] = expression_values
    return cell_names, exp_by_gene

## Save to firestore
## Abstract this out
def add_data_to_firestore(exp_by_gene, cells, matrix_path):
    batch_start_time = time.time()
    logger.info('Start batch insert')
    batch = db.batch()
    for gene, exp_vals in exp_by_gene.items():
        doc_ref = db.collection('genes').document(gene)
        doc_data = {
            'source_file_name': matrix_path,
            'data_arrays': {
                'cell_names': cells,
                'expression_scores': exp_vals # TODO: Align names: "score" or "value"?
            }
        }
        batch.set(doc_ref, doc_data)
    batch.commit()
    logger.info('End batch insert')
    batch_time = str(round(time.time() - batch_start_time))
    logger.info('Batch insert time: ' + batch_time + ' s')
    time.sleep(2)

def chunk(data, SIZE=5):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}

logger = get_logger('', 'ingest')
cells, exp_by_gene = map_genes_to_expression_values(args.matrix_path)
# for genes in chunk(exp_by_gene):
add_data_to_firestore(exp_by_gene, cells, args.matrix_path)


logger.info('End dense_ingest.py')
total_time = str(round(time.time() - start_time))
logger.info('Time: ' + total_time + ' s')