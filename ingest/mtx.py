"""
Module for ingesting MTX files

DESCRIPTION
This module provides extract and transforms function for a MTX file bundle. A
MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

"""

import datetime
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

from bson.objectid import ObjectId

try:
    from expression_files import GeneExpression
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in
    # single_cell_portal code
    from .expression_files import GeneExpression
    from .ingest_files import IngestFiles


class MTXIngestor(GeneExpression):
    ALLOWED_FILE_TYPES = ["text/tab-separated-values"]

    def matches_file_type(file_type):
        return "mtx" == file_type

    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
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
        self.mtx_file = mtx_ingest_file.resolve_path(mtx_path)[0]
        # Only known way to traverse through zipped files w/o unzipping
        next(self.mtx_file)
        next(self.mtx_file)
        # A list ['N', 'K', 'M'] that represents a gene-barcode matrix where N
        # is the gene index, M is the barcode index, and K is the expression
        # score for the given gene index
        self.mtx_description = next(self.mtx_file)

    @staticmethod
    def is_sorted(idx: int, visited_expression_idx: int):
        last_visited_idx = visited_expression_idx[-1]
        if idx not in visited_expression_idx:
            if idx == (last_visited_idx + 1):
                return True
            else:
                return False
        else:
            if idx == last_visited_idx:
                return True
            else:
                return False

    def create_data_array(
        self, gene: str, array_type: str, values: List[str], linear_data_type: str, id
    ) -> List:
        data_arrays = []
        for model in self.set_data_array(
            gene, array_type, values, linear_data_type, id
        ):
            data_arrays.append(model)
        return data_arrays

    def transform(self):
        start_time = datetime.datetime.now()
        num_processed = 0
        last_idx = 0
        gene_models = []
        data_arrays = []
        exp_cells = []
        exp_scores = []
        visited_expression_idx = [0]
        model_id = None

        # All cells that were observed in expression matrix
        data_arrays.extend(
            self.create_data_array(
                f"{self.cluster_name} Cells",
                "cells",
                self.cells,
                "Study",
                self.study_file_id,
            )
        )
        for row in self.mtx_file:
            raw_gene_idx, raw_barcode_idx, raw_exp_score = row.split()
            current_idx = int(raw_gene_idx)
            gene_id, gene = self.genes[current_idx - 1].split("\t")
            if current_idx != last_idx:
                if MTXIngestor.is_sorted(current_idx, visited_expression_idx):
                    visited_expression_idx.append(current_idx)
                    # Create data arrays from prior gene
                    if last_idx != 0:
                        # Cell names that had significant (i.e. non-zero)
                        # expression for gene
                        dr_models = self.create_data_array(
                            gene, f"{gene} Cells", exp_cells, "Study", model_id
                        )
                        data_arrays.extend(dr_models)
                        # Significant (i.e. non-zero) expression values for
                        # gene
                        dr_models = self.create_data_array(
                            gene, f"{gene} Expression", exp_scores, "Gene", model_id
                        )
                        data_arrays.extend(dr_models)
                        # Reset variables so values will be associated w/new
                        # gene
                        exp_cells = []
                        exp_scores = []
                    if len(data_arrays) > 1_000:
                        num_processed += len(gene_models)
                        print(
                            f"Processed {num_processed} models, "
                            f"{str(datetime.datetime.now() - start_time)} "
                            f"elapsed"
                        )
                        yield gene_models, data_arrays
                        gene_models = []
                        data_arrays = []
                    model_id = ObjectId()
                    # Add current gene's gene model
                    gene_model = self.create_gene_model(
                        gene, self.study_file_id, self.study_id, gene_id, model_id
                    )
                    gene_models.append(gene_model)
                    last_idx = current_idx
                else:
                    raise ValueError("MTX file must be sorted")
            exp_cell = self.cells[int(raw_barcode_idx) - 1]
            exp_score = round(float(raw_exp_score), 3)
            exp_cells.append(exp_cell)
            exp_scores.append(exp_score)
        yield gene_models, data_arrays
        num_processed += len(gene_models)
        print(
            f"Processed {num_processed} models, "
            f"{str(datetime.datetime.now() - start_time)} "
            f"elapsed"
        )

    def execute_ingest(self):
        self.extract_feature_barcode_matrices()
        for gene, data_arrays in self.transform():
            yield gene, data_arrays

    def extract_feature_barcode_matrices(self):
        """
        Sets relevant iterables for the gene and barcode file of the MTX bundle
        """
        self.genes = [g.strip().strip('"') for g in self.genes_file.readlines()]
        self.cells = [c.strip().strip('"') for c in self.barcodes_file.readlines()]
