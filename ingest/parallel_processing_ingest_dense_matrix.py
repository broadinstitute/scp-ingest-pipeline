# """THIS IS A WORK IN PROGRESS
# Command-line interface for ingesting dense matrixes into firestore with paraelle processing

# DESCRIPTION
# This CLI uses paraelle processing to map expression data from a dense matrix.
# Genes are mapped to expression values and cell  names and inputted into firestore

# PREREQUISITES
# You must have google could firestore installed, authenticated
#  configured. !!!!!! Also make sure to run command:
#  export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES !!!!!!!!

# EXAMPLES
# # Takes dense file and stores it into firestore
# python parallel_processing_ingest_dense_matrix.py <file path>
# $ python parallel_processing_ingest_dense_matrix.py --file-path ../tests/data/dense_matrix_19_genes_100k_cells.txt

# """
# import argparse
# import multiprocessing
# import sys
# import time
# from itertools import islice
# from signal import SIG_DFL, SIGPIPE, signal

# from google.cloud import firestore

# # Broken pipe error is happening because main program finished before processes causing error
# # BrokenPipeError: [Errno 32] Broken pipe. This is a work around for now
# signal(SIGPIPE, SIG_DFL)

# db = firestore.Client()


# parser = argparse.ArgumentParser(
#     prog='parallel_processing_ingest_dense_matrix_.py',
#     description=__doc__,
#     formatter_class=argparse.RawDescriptionHelpFormatter,
# )

# # Positional argument
# parser.add_argument(
#     "--file-path", required=True, help='Absolute or relative path to dense matrix'
# )

# args = parser.parse_args()


# def determine_size_of_chunks(file_object, file_path):
#     """Determines amount of lines to be returned per chunk
#     Ensures field value, in this case a line, is less than 1,048,487 bytes and,
#     eventually, API request size is less than 10 MiB. If field value is more
#     than 1,048,487 bytes, logic will be needed to figure out how chunk per line
#     basis.

#     TODO:
#     Abstract out function. Needs to handle different ways to chunk based on file type and size restrants
#     """

#     line = file_object.readline()

#     # file_size will be eventually used to track loading and transforming progress
#     # file_size = os.path.getsize(file_path)
#     size_of_single_row = sys.getsizeof(line)
#     print("Size of one line is %i " % size_of_single_row)
#     amount_of_lines = 1048487 / size_of_single_row
#     print("Amount of rows in 1 MiB %s" % amount_of_lines)
#     if amount_of_lines > 500:
#         return 500
#     else:
#         return int(amount_of_lines)


# def extract(file, number_of_lines):
#     # skip first line for dense matrix
#     file.readline()
#     while True:
#         next_lines = list(islice(file, number_of_lines))
#         if not next_lines:
#             # Let's worker function that there are no more extracted data
#             return []
#             break
#         yield next_lines


# # worker function


# def transform(process_name, extracted_data, transformed_data):
#     print('[%s] transform routine starts' % process_name)
#     gene_list = []
#     model = {}

#     while True:
#         try:
#             new_value, expression_attr = extracted_data.get()

#             if not new_value:
#                 print('[%s] transform routine quits' % process_name)
#                 # Indicate transformation function is finished and signal load worker function to quit
#                 transformed_data.put([])
#                 break
#             else:
#                 print('[%s] is transforming data' % process_name)
#                 # Transform data into datamodel
#                 for line in new_value:
#                     compute = line.rstrip('\n').split(',')
#                     model[compute[0]] = {}
#                     model[compute[0]]['gene_expression'] = {}
#                     model[compute[0]]['gene_expression']['cell_names'] = expression_attr
#                     # change this line to alter amount of scores. Right now first
#                     # 1k scores and genes are loaded to data model because  we
#                     # cannot yet ingest matrices with many cells,
#                     # due to one of the Firestore constraints
#                     scores = [float(x) for x in compute[1:1000]]
#                     model[compute[0]]['gene_expression']['expression_scores'] = scores
#                     gene_list.append(model)
#                 transformed_data.put(model)
#             extracted_data.task_done()
#         except EOFError as error:
#             # signal transform worker to quit if the queue is empty
#             print('[%s] transform routine quits' % process_name)
#             print(error)
#             break

#     return


# def load(process_name, transformed_data):
#     print('[%s] load routine starts' % process_name)

#     while True:
#         try:
#             list_of_transformed_data = transformed_data.get()

#             if not list_of_transformed_data:
#                 print('Process Load [%s]  routine quits' % process_name)
#                 break

#             else:
#                 batch = db.batch()
#                 for key, val in list_of_transformed_data.items():
#                     print(' loading %s ' % key)
#                     doc_ref = db.collection("gene").document(key)
#                     batch.set(doc_ref, val)
#                 batch.commit()
#                 time.sleep(0.2)
#             transformed_data.task_done()
#         except EOFError as error:
#             # signal load worker function to quit if the queue is empty
#             print('Process Load [%s]  routine quits' % process_name)
#             print(error)
#             break

#     return


# if __name__ == "__main__":
#     # Define manager
#     manager = multiprocessing.Manager()

#     # Define a queue for extracted  and transformed data
#     extracted_data = manager.Queue()
#     transformed_data = manager.Queue()

#     # Decide amount of processes to run. There are 2 processes so amount is divided by 2
#     # 1 is subtracted from count for safe measures
#     num_processes = int((multiprocessing.cpu_count() - 1) / 2)
#     processes = []
#     expression_attr = []

#     # Initiate the worker processes
#     for i in range(num_processes):
#         # Set process name
#         transform_process_name = 'Transform P%i' % i
#         load_process_name = 'Load P%i' % i

#         # Create the process, and connect it to the worker function
#         transform_process = multiprocessing.Process(
#             target=transform,
#             args=(transform_process_name, extracted_data, transformed_data),
#         )
#         load_process = multiprocessing.Process(
#             target=load, args=(load_process_name, transformed_data)
#         )

#         # Start the process
#         transform_process.start()
#         load_process.start()

#     # Fill extract queue
#     with open(args.file_path, 'r') as fname:
#         # change this line to alter amount of cells. Right now first 1k are loaded
#         # because  we cannot yet ingest matrices with many cells,
#         # due to one of the Firestore constraints
#         cell_names = fname.readline().rstrip().split(',')[1:1000]
#         print("Size of cell names is %i " % sys.getsizeof(cell_names))
#         size_of_chunks = determine_size_of_chunks(fname, args.file_path)
#         fname.seek(0, 0)
#         for data in extract(fname, size_of_chunks):
#             extracted_data.put((data, cell_names))
#             time.sleep(0.1)
