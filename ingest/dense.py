"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""
import argparse
import os
import sys
from itertools import islice
from typing import *

import numpy as np
from gene_data_model import Gene
from ingest_files import IngestFiles


class Dense():
    def __init__(self, file_path):
<<<<<<< Updated upstream
        self.file = open(file_path, 'r')
        self.cell_names = self.file.readline().replace('"', '').split(',')[1:]

        self.file_name = file_path.strip(".")

    def extract(self, size: int = 500) -> List[str]:
        """Extracts lines from dense matrix.

        Args:
            size : int
                The amount of lines returned per chunk

        Returns:
                next_lines : List[str]
                    A list (chunk) of rows from a dense matrix.
        """
        while True:
            next_lines = list(islice(self.file, size))
            if not next_lines:
                break
            yield next_lines
=======
        self.ALLOWED_FILE_TYPES = ['text/csv',
                                   'text/plain', 'text/tab-separated-values']
        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.cell_names = self.get_next_line(increase_line_count=False)[1:]
>>>>>>> Stashed changes

    def transform_expression_data_by_gene(self, expression_scores: List[str]) -> Gene:
        """Transforms dense matrix into firestore data model for genes.

        Args:
            lines : List[str]
                Lines from dense matrix file

        Returns:
                transformed_data : List[Gene]
                A list of Gene objects
        """
<<<<<<< Updated upstream
        transformed_data = []
        print(f'starting ingesting {len(lines)} lines')
        for line in lines:
            expression_scores = line.rstrip('\n').split(',')
            compute = line.rstrip('\n').split(',')
            gene_model = Gene(compute[0], source_file_type="Dense",
                              expression_scores=expression_scores[1:],
                              cell_names=self.cell_names)
            transformed_data.append(gene_model)
        return transformed_data
=======
        gene_model = Gene(expression_scores[0], source_file_type="Dense",
                          expression_scores=expression_scores[1:],
                          cell_names=self.cell_names)
        return gene_model
>>>>>>> Stashed changes

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.file.close()
