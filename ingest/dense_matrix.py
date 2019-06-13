"""Command-line interface for ingesting dense matrixes into firestore

DESCRIPTION
Genes are mapped to expression values and cell  names and inputted into firestore

PREREQUISITES
Must have python 3.6 or higher.

EXAMPLES
# Takes dense file and stores it into firestore
$ python dense_matrix.py <file path>
$ python dense_matrix.py ../tests/data/dense_matrix_19_genes_100k_cells.txt
"""
import argparse
from typing import *
import numpy as np
import sys
import os
from itertools import islice

from gene_data_model import Gene
import ingest

class Dense():
    def __init__(self, file_path):
        if not os.path.exists(file_path):
		          raise IOError(f"File '{file_path}' not found")
        self.file = open(file_path,'r')
        self.cell_names = self.file.readline()[1:1000]
        self.file_name, self.filetype = os.path.splitext(file_path)

    def extract(self, size = 500):
        while True:
            next_lines = list(islice(self.file, size))
            if not next_lines:
                break
            yield next_lines

    def transform(self, *lines):
        """Transforms dense matrix into firestore data model

        Args:
            lines : List[str]
                Lines from dense matrix file

        Returns:
        ------
		    transformed_data : List[Gene]
                A list of Gene objects
        """
        transformed_data = []
        for line in lines:
            compute = line.rstrip('\n').split(',')
            expression_scores = [float(x) for x in compute[1:1000]]
            gene_model = Gene(compute[0], self.file_name, self.filetype,
            expression_scores = expression_scores,
            cell_names = self.cell_names)
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
        prog = 'dense_matrix.py',
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )

    #Positional argument
    parser.add_argument(
        "dense_matrix",
        help ='Absolute or relative path to expression file'
    )
    args = parser.parse_args()

    dense_matrix_object = Dense(args.dense_matrix)
    dense_matrix_object.ingest()
    dense_matrix_object.loom.close()
