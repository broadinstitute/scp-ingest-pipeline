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

try:
    from expression_files import GeneExpression
    from monitor import setup_logger, log
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .expression_files import GeneExpression


class Mtx(GeneExpression):
    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)
        self.genes_file = open(kwargs.pop("gene_file"))
        self.barcodes_file = open(kwargs.pop("barcode_file"))
        self.mtx_path = mtx_path
        self.matrix_params = kwargs
        self.exp_by_gene = {}

    def extract(self):
        """Sets relevant iterables for each file of the MTX bundle
        """
        self.matrix_file = scipy.io.mmread(self.mtx_path)
        self.genes = [g.strip() for g in self.genes_file.readlines()]
        self.cells = [c.strip() for c in self.barcodes_file.readlines()]

    def transform(self):
        """Transforms matrix gene data model
        """
        # Named tuple that has orignial string format of a gene's name, 'gene_name'
        # and it's corresponding model 'gene_model'
        GeneModelData = collections.namedtuple('GeneModel', ['gene_name', 'gene_model'])
        # Named tuple to hold expression data for a single gene
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
                # Append new score to 'expression_scores' key in GeneModelData
                self.exp_by_gene[gene].expression_scores.append(exp_score)
                # Append new cell name to 'cell_names' key in GeneModelData
                self.exp_by_gene[gene].cell_names.append(cell_name)
            else:
                self.exp_by_gene[gene] = GeneExpressionValues([cell_name], [exp_score])
                self.info_logger.info(
                    f'Creating model for {gene} ', extra=self.extra_log_params
                )
                yield GeneModelData(
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
        self,
        linear_data_id,
        unformatted_gene_name,
        gene_name,
        create_cell_DataArray=False,
    ):
        if create_cell_DataArray:
            self.info_logger.info(
                f'Creating cell data array for gene :{unformatted_gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_cells(self.cells, linear_data_id)
        else:
            self.info_logger.info(
                f'Creating cell names data array for gene:{unformatted_gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_gene_cell_names(
                unformatted_gene_name,
                linear_data_id,
                self.exp_by_gene[unformatted_gene_name].cell_names,
            )
            self.info_logger.info(
                f'Creating gene expression data array for gene:{unformatted_gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_gene_expression_values(
                unformatted_gene_name,
                linear_data_id,
                self.exp_by_gene[unformatted_gene_name].expression_scores,
            )

    def close(self):
        """Closes file
        """
        self.genes_file.close()
        self.barcodes_file.close()
