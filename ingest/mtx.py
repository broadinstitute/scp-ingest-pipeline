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
from bson.objectid import ObjectId
import copy

try:
    from expression_files import GeneExpression
    from ingest_files import IngestFiles
    from monitor import trace
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .expression_files import GeneExpression
    from .ingest_files import IngestFiles
    from .monitor import trace


class Mtx(GeneExpression):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        self.tracer = kwargs.pop("tracer")
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)

        genes_path = kwargs.pop("gene_file")
        genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
        self.genes_file, genes_local_path = genes_ingest_file.resolve_path(genes_path)

        barcodes_path = kwargs.pop("barcode_file")
        barcodes_ingest_file = IngestFiles(barcodes_path, self.ALLOWED_FILE_TYPES)
        self.barcodes_file, barcodes_local_path = barcodes_ingest_file.resolve_path(
            barcodes_path
        )

        mtx_ingest_file = IngestFiles(mtx_path, self.ALLOWED_FILE_TYPES)
        mtx_file, self.mtx_local_path = mtx_ingest_file.resolve_path(mtx_path)

        self.matrix_params = kwargs
        self.exp_by_gene = {}

    @trace
    def extract(self):
        """Sets relevant iterables for each file of the MTX bundle
        """
        self.matrix_file = scipy.io.mmread(self.mtx_local_path)
        self.genes = [g.strip() for g in self.genes_file.readlines()]
        self.cells = [c.strip() for c in self.barcodes_file.readlines()]

    @trace
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
            if gene in self.exp_by_gene.keys():
                # Append new score to 'expression_scores' key in GeneModelData
                self.exp_by_gene[gene].expression_scores.append(exp_score)
                # Append new cell name to 'cell_names' key in GeneModelData
                self.exp_by_gene[gene].cell_names.append(cell_name)
            else:
                self.exp_by_gene[gene] = copy.copy(
                    GeneExpressionValues([exp_score], [cell_name])
                )
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
                            '_id': ObjectId(),
                        }
                    ),
                )

    @trace
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
