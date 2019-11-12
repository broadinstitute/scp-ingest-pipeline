"""
Module for ingesting MTX files

DESCRIPTION
This module provides extract and transforms function for gene expression data for
an MTX file bundle. An MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

PREREQUISITES
Must have python 3.6 or higher.
"""

from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
import collections
import scipy.io
from expression_files import GeneExpression


class Mtx(GeneExpression):
    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)
        self.genes_file = open(kwargs.pop("gene_file"))
        self.barcodes_file = open(kwargs.pop("barcode_file"))
        self.mtx_path = mtx_path
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.matrix_params = kwargs
        self.exp_by_gene = {}

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

    def transform(self):
        """Transforms dense matrix into firestore data model for genes.

        Args:
            None

        Returns:
            transformed_data : List[Gene]
                A list of Gene objects
        """
        GeneModel = collections.namedtuple('GeneModel', ['gene_name', 'gene_model'])

        GeneExpressionValues = collections.namedtuple(
            'GeneExpressionValues', ['expression_scores', 'cell_names']
        )
        for raw_gene_idx, raw_barcode_idx, raw_exp_score in zip(
            self.matrix_file.row, self.matrix_file.col, self.matrix_file.data
        ):
            gene_id, gene = self.genes[int(raw_gene_idx)].split('\t')
            cell_name = self.cells[int(raw_barcode_idx)]
            exp_score = round(float(raw_exp_score), 3)
            if gene in self.exp_by_gene:
                # Append new score to 'expression_scores' key in Gene object
                self.exp_by_gene[gene].expression_scores.append(exp_score)
                self.exp_by_gene[gene].cell_names.append(cell_name)
            else:
                self.exp_by_gene[gene] = GeneExpressionValues([cell_name], [exp_score])
                yield GeneModel(
                    gene,
                    self.Model(
                        {
                            'name': gene,
                            'searchable_name': gene.lower(),
                            'study_file_id': self.study_file_id,
                            'study_id': self.study_id,
                            'gene_id': gene_id,
                        }
                    ),
                )

    def set_data_array(
        self, unformatted_gene_name, name, linear_data_id, create_cell_DataArray=False
    ):
        print(self.exp_by_gene[unformatted_gene_name])
        if create_cell_DataArray:
            yield self.set_data_array_cells(self.cells, linear_data_id)
        else:
            yield self.set_data_array_gene_cell_names(
                name, linear_data_id, self.exp_by_gene[unformatted_gene_name].cell_names
            )
            yield self.set_data_array_gene_expression_values(
                name,
                linear_data_id,
                self.exp_by_gene[unformatted_gene_name].expression_scores,
            )

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.genes_file.close()
        self.barcodes_file.close()
