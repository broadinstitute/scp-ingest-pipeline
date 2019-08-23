"""Module for ingesting MTX

DESCRIPTION
This module provides extract and transforms function for gene expression data for
an MTX file bundle. An MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

PREREQUISITES
Must have python 3.6 or higher.
"""

from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

import scipy.io
from gene_data_model import Gene


class Mtx:
    def __init__(self, mtx_path: str, genes_path: str, barcodes_path: str):
        self.genes_file = open(genes_path)
        self.barcodes_file = open(barcodes_path)
        self.mtx_path = mtx_path

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
        for raw_gene_idx, raw_barcode_idx, raw_exp_score in zip(
            self.matrix_file.row, self.matrix_file.col, self.matrix_file.data
        ):
            gene_id, gene = self.genes[int(raw_gene_idx)].split('\t')
            cell_name = self.cells[int(raw_barcode_idx)]
            exp_score = round(float(raw_exp_score), 3)
            if gene in exp_by_gene:
                # Append new score to 'expression_scores' key in Gene object
                exp_by_gene[gene].expression_scores.append(exp_score)
                exp_by_gene[gene].cell_names.append(cell_name)
            else:
                # Create new key value pair with value being Gene object
                exp_by_gene[gene] = Gene(
                    gene,
                    'Mtx',
                    gene_id=gene_id,
                    cell_names=[cell_name],
                    expression_scores=[exp_score],
                    check_for_zero_values=False,
                )
        return exp_by_gene.values()

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.genes_file.close()
        self.barcodes_file.close()
