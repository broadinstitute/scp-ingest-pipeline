"""Command-line interface for ingesting dense matrixes into firestore

DESCRIPTION
Genes are mapped to expression values and cell  names and inputted into firestore

PREREQUISITES
You must have google could firestore installed, authenticated
 configured. 

EXAMPLES
# Takes dense file and stores it into firestore
python parallel_processing_ingest_dense_matrix.py <file path>
$ python parallel_processing_ingest_dense_matrix.py --file-path ../tests/data/dense_matrix_19_genes_100k_cells.txt
"""
import argparse
from typing import *
import numpy as np
import sys
import os

from gene_data_model import Gene
import ingest

class Dense():
    def __init__(self, file_path):
        if not os.path.exists(file_path):
		          raise IOError(f"File '{file_path}' not found")
        self.file = open(file_path,'r')
        self.cell_names = self.file.readline()[1:1000]
        self.file_name, self.filetype = os.path.splitext(file_path)

    def extract(file, size=500):
        while True:
            next_lines = list(islice(file, number_of_lines))
            if not next_lines:
                break
            yield next_lines

    def transform(lines):
        transformed_data = []
        for line in lines:
            compute = line.rstrip('\n').split(',')
            expression_scores = [float(x) for x in compute[1:1000]]
            gene_model = Gene(compute[0], self.file_name, self.filetype,
            expression_scores = expression_scores,
            cell_names = self.cell_name)
            transformed_data.append(gene_model.gene)
        return transformed_data

    def ingest(self) -> None:
        """ Ingests dense matrix file via Ingest_Service

        Args:
            Nothing
        Returns:
        ------
            Nothing
        """
        ingest_service = ingest.connect(self.extract, self.transform)
        ingest_service.ingest()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog ='dense_matrix.py',
        description= __doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    #Positional argument
    parser.add_argument(
        "dense-matrix",
        help='Absolute or relative path to expression file'
    )
    args = parser.parse_args()

    dense_matrix_object = Dense(args.dense_matrix)
    dense_matrix_object.ingest()
    dense_matrix_object.loom.close()
