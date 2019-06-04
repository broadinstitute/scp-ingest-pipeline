"""Command-line interface for ingesting dense matrixes into firestore with paraelle processing

DESCRIPTION
This CLI uses paraelle processing to map expression data from a dense matrix.
Genes are mapped to expression values and cell  names and inputted into firestore

PREREQUISITES
You must have google could firestore installed, authenticated
 configured. Also make sure to run command:
 export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

EXAMPLES
# Takes dense file and stores it into firestore
$ python paraelle_processing_dense_matrix.py  200k_subsample_4k_PMBC.loom
"""
import multiprocessing
import argparse
import time
import sys
import os
from itertools import islice
from sys import getsizeof
import linecache
import queue

from google.cloud import firestore
import loompy
import numpy as np

db = firestore.Client()

# define line_number
line_number = 3

parser = argparse.ArgumentParser(
    prog ='parallel_processing_expressionData.py'
)

parser.add_argument(
    "file_name", default=None,
    help='Path to file'
)

args = parser.parse_args()
read_shutdown_value = False
write_shutdown_value = False


def determine_size_of_chunks(file):
    line = linecache.getline(args.file_name, line_number).split()
    file_size = os.path.getsize(args.file_name)
    size_of_single_row = sys.getsizeof(line)*.000001
    amount_of_lines = 1/size_of_single_row
    print("Amount of rows in 1MB %s" % amount_of_lines)
    print(amount_of_lines)
    if amount_of_lines > 500:
        return 500
    else:
        return int(amount_of_lines)

def extract(file, number_of_lines):
    it = iter(file)
    while True:
        next_lines=list(islice(file, 5))
        yield next_lines
        if not next_lines:
            # Let's worker function that there are no more extracted data
            return []
            break

##worker function
def transform(process_name, extracted_data, transformed_data,expression_attr):
    print('[%s] transform routine starts' % process_name)
    gene_list = []
    model = dict()

    while True:
        try:
            new_value = extracted_data.get()

            if len(new_value) < 1:
                print('[%s] transform routine quits' % process_name)
                # Indicate finished and signal load worker function to quit
                transformed_data.put([])
                break
            else:
                #Transform data into datamodel
                for line in new_value:
                    compute=line.rstrip('\n').split(',')
                    model[compute[0]]= {}
                    model[compute[0]]['expression_values'] = compute[1:2000]
                    model[compute[0]]['cell_names'] = expression_attr
                    gene_list.append(model)
                    print('[%s] Added values: %s' % (process_name, compute[0]))
                transformed_data.put(model)
                time.sleep(.1)

        except EOFError as error:
            # Output expected EOFErrors.
            print('[%s] transform routine quits' % process_name)
            break

    return


##Abstract this out
def load(process_name, transformed_data):
    print('[%s] load routine starts' % process_name)

    while True:
        try:
            list_of_transformed_data = transformed_data.get()

            if len(list_of_transformed_data)< 1:
                print('Process Load [%s]  routine quits' % process_name)
                break

            else:
                batch = db.batch()
                for key,val in list_of_transformed_data.items():
                    doc_ref = db.collection("dense_matrix_gene").document(key)
                    batch.set(doc_ref, val)
                batch.commit()
                time.sleep(.2)

        except EOFError as error:
            # Output expected EOFErrors.
            print('Process Load [%s]  routine quits' % process_name)
            break
    return


if __name__ == "__main__":
    # Define IPC manager
    manager = multiprocessing.Manager()

    # Define a list (queue) for extraced data and transformed data
    extracted_data = manager.Queue()
    transformed_data = manager.Queue()

    num_processes= int((multiprocessing.cpu_count()-1)/2)
    pool = multiprocessing.Pool(processes=num_processes)
    processes = []
    expression_attr = []

    # Initiate the worker processes
    for i in range(num_processes):
        # Set process name
        transform_process_name = 'Transform P%i' % i
        load_process_name = 'Load P%i' % i

        # Create the process, and connect it to the worker function
        transform_process = multiprocessing.Process(target=transform, args=(transform_process_name,extracted_data,transformed_data, expression_attr ))
        load_process = multiprocessing.Process(target=load, args=(load_process_name, transformed_data ))

        # Start the process
        transform_process.start()
        load_process.start()

    # Fill extract queue
    with open(args.file_name,'r') as fname:
        expression_attr = fname.readline().rstrip().split(',')[1:]
        size_of_chunks = determine_size_of_chunks(fname)
        for data in extract(fname, size_of_chunks):
            extracted_data.put(data)
    # Wait while the workers process
    time.sleep(5)
    pool.close()
    pool.join()
    transform_process.join()
