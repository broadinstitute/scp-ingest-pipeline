"""
Module for ingesting MTX files

DESCRIPTION
This module provides extract and transforms function for a MTX file bundle. A
MTX file bundle consists of A) an .mtx file in Matrix
Market matrix coordinate format, B) a genes.tsv file, and a barcodes.tsv file.
These are commonly provided from 10x Genomics v2.

"""

import datetime
from typing import Dict, Generator, List, Tuple  # noqa: F401
import sys

from bson.objectid import ObjectId

try:
    from expression_files import GeneExpression

    sys.path.append("../ingest")
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in
    # single_cell_portal code
    from .expression_files import GeneExpression

    sys.path.append("../ingest")
    from ..ingest_files import IngestFiles


class MTXIngestor(GeneExpression):
    ALLOWED_FILE_TYPES = ["text/tab-separated-values"]

    @staticmethod
    def matches_file_type(file_type):
        return "mtx" == file_type

    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)

        mtx_ingest_file = IngestFiles(mtx_path, self.ALLOWED_FILE_TYPES)
        self.mtx_file = mtx_ingest_file.resolve_path(mtx_path)[0]

        genes_path = kwargs.pop("gene_file")
        genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
        self.genes_file = genes_ingest_file.resolve_path(genes_path)[0]

        barcodes_path = kwargs.pop("barcode_file")
        barcodes_ingest_file = IngestFiles(barcodes_path, self.ALLOWED_FILE_TYPES)
        self.barcodes_file = barcodes_ingest_file.resolve_path(barcodes_path)[0]

        # A list ['N', 'K', 'M'] that represents a gene-barcode matrix where N
        # is the gene index, M is the barcode index, and K is the expression
        # score for the given gene index
        self.mtx_dimensions: List[int] = MTXIngestor.get_mtx_dimensions(self.mtx_file)

    @staticmethod
    def check_valid(cell_names: List[str], barcodes, genes):
        error_messages = []

        try:
            MTXIngestor.check_bundle(barcodes, genes, cell_names)
        except ValueError as v:
            error_messages.append(str(v))
        try:
            MTXIngestor.check_duplicate_genes(genes)
        except ValueError as v:
            error_messages.append(str(v))
        try:
            GeneExpression.check_unique_cells(cell_names)
        except ValueError as v:
            error_messages.append(str(v))

        if len(error_messages) > 0:
            raise ValueError("; ".join(error_messages))
        return True

    @staticmethod
    def check_bundle(barcodes, genes, mtx_dimensions):
        """Confirms barcode and gene files have expected length of values"""
        expected_genes = mtx_dimensions[0]
        expected_barcodes = mtx_dimensions[1]
        actual_genes = len(barcodes)
        actual_barcodes = len(genes)

        if actual_barcodes == expected_barcodes and actual_genes == expected_genes:
            msg = (
                f"Expected {expected_barcodes} cells and {expected_genes}. "
                f"Got {actual_barcodes} cells and {actual_genes}"
            )
            raise ValueError(msg)
        return True

    @staticmethod
    def check_duplicate_genes(gene_names: List):
        unique_gene_names: List[str] = set(gene_names)
        if len(unique_gene_names) != len(gene_names):
            amount_of_duplicates = unique_gene_names - len(gene_names)
            msg = (
                "Duplicate header values are not allowed."
                f"There are {amount_of_duplicates} duplicates"
            )
            raise ValueError(msg)
        return True

    @staticmethod
    def get_mtx_dimensions(file_handler) -> List:
        for line in file_handler:
            if not line.startswith("%"):
                mtx_dimensions: List[str] = line.strip().split()
                # Convert values in mtx_dimensions to int
                return list(map(int, mtx_dimensions))

    @staticmethod
    def is_sorted(idx: int, visited_expression_idx: List[int]):
        last_visited_idx = visited_expression_idx[-1]
        if idx not in visited_expression_idx:
            if idx == (last_visited_idx + 1):
                return True
            else:
                raise ValueError("MTX file must be sorted")
        else:
            if idx == last_visited_idx:
                return True
            else:
                raise ValueError("MTX file must be sorted")

    def create_data_array(
        self, gene: str, array_type: str, values: List[str], linear_data_type: str, id
    ) -> List:
        data_arrays = []
        for model in GeneExpression.set_data_array(
            gene,
            self.cluster_name,
            array_type,
            values,
            linear_data_type,
            id,
            self.study_id,
            self.study_file_id,
        ):
            data_arrays.append(model)
        return data_arrays

    def execute_ingest(self):
        self.extract_feature_barcode_matrices()
        if MTXIngestor.check_valid(self.cells, self.genes, self.mtx_dimensions):
            for gene_docs, data_array_documents in self.transform:
                self.load(gene_docs, data_array_documents)

    def extract_feature_barcode_matrices(self):
        """
        Sets relevant iterables for the gene and barcode file of the MTX bundle
        """
        self.genes: List[str] = [
            g.strip().strip('"') for g in self.genes_file.readlines()
        ]
        self.cells: List[str] = [
            c.strip().strip('"') for c in self.barcodes_file.readlines()
        ]

    @property
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

        # All observed cells
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
                is_sorted = MTXIngestor.is_sorted(current_idx, visited_expression_idx)
                if not is_sorted:
                    raise ValueError("MTX file must be sorted")
                else:
                    visited_expression_idx.append(current_idx)
                    # Create data arrays from prior gene
                    if last_idx != 0:
                        # Data array for cell names
                        da_cells = self.create_data_array(
                            gene, f"{gene} Cells", exp_cells, "Study", model_id
                        )
                        data_arrays.extend(da_cells)
                        # Data array for expression values
                        da_exp = self.create_data_array(
                            gene, f"{gene} Expression", exp_scores, "Gene", model_id
                        )
                        data_arrays.extend(da_exp)
                        # Reset variables so values will be associated w/new
                        # gene
                        exp_cells = []
                        exp_scores = []
                if len(data_arrays) > 1_000:
                    yield gene_models, data_arrays
                    num_processed += len(gene_models)
                    print(
                        f"Processed {num_processed} genes. "
                        f"{str(datetime.datetime.now() - start_time)} "
                        f"elapsed"
                    )
                    num_processed += len(gene_models)
                    gene_models = []
                    data_arrays = []
                model_id = ObjectId()
                # Add current gene's gene model
                gene_model = GeneExpression.create_gene_model(
                    gene, self.study_file_id, self.study_id, gene_id, model_id
                )
                gene_models.append(gene_model)
                last_idx = current_idx
            exp_cell = self.cells[int(raw_barcode_idx) - 1]
            exp_score = round(float(raw_exp_score), 3)
            exp_cells.append(exp_cell)
            exp_scores.append(exp_score)
        if gene_models:
            yield gene_models, data_arrays
            num_processed += len(gene_models)
            print(
                f"Processed {num_processed} genes. "
                f"{str(datetime.datetime.now() - start_time)} "
                f"elapsed"
            )
