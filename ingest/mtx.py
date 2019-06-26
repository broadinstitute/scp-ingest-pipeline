"""Module for ingesting MTX

DESCRIPTION
This module provides extract and transforms function for gene expression data for
an MTX file bundle. An MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

PREREQUISITES
Must have python 3.6 or higher.
"""
import os
import sys
from typing import Dict, Generator, List, Tuple, Union

import scipy.io
from gene_data_model import Gene


class Mtx:
    def __init__(self, mtx_path, bundle_paths):
        genes_path, barcodes_path = bundle_paths
        self.genes_file = open(genes_path)
        self.barcodes_file = open(barcodes_path)
        self.mtx_path, self.source_file_type = os.path.splitext(
            mtx_path)

    def extract(self):
        """Sets relevant iterables for each file of the MTX bundle

        Args:
            None

        Returns:
            None
        """
        self.matrix_file = scipy.io.mmread(self.mtx_path)
        self.genes = [g.strip() for g in self.genes_file.readlines()]
        self.cells = [c.strip() for c in self.barcodes_file.readlines()]

    def transform_expression_data_by_gene(self) -> List[Gene]:
        """Transforms dense matrix into firestore data model for genes.

        Args:
            None

        Returns:
            transformed_data : List[Gene]
                A list of Gene objects
        """
        exp_by_gene = {}
        for i, j, exp_score in zip(self.matrix_file.row, self.matrix_file.col, self.matrix_file.data):
            gene_id, gene = self.genes[i].split('\t')
            exp_score = float(exp_score)
            if gene in exp_by_gene:
                # Append new score to 'expression_scores' key in Gene object
                exp_by_gene[gene].get_expression_scores().append(exp_score)
            else:
                # Create new key value pair with value being Gene object
                exp_by_gene[gene] = Gene(gene, self.mtx_path,
                                         self.source_file_type, gene_id=gene_id,
                                         cell_names=self.cells,
                                         expression_scores=[exp_score])
        # Get 'gene' value from Gene Object
        transformed_data = {k: v.gene for k, v in exp_by_gene.items()}
        return transformed_data.values()

    def close(self):
        self.genes_file.close()
        self.barcodes_file.close()
