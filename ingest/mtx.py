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

import collections
import copy
import datetime
import linecache
import os
import subprocess
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

import scipy.io
from bson.objectid import ObjectId

try:
    from expression_files import GeneExpression
    from ingest_files import IngestFiles
    from monitor import trace
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .expression_files import GeneExpression
    from .ingest_files import IngestFiles
    from .monitor import trace


class MTXIngestor(GeneExpression):
    ALLOWED_FILE_TYPES = ["text/tab-separated-values"]

    def matches_file_type(file_type):
        return 'mtx' == file_type

    def __init__(self, mtx_path: str, study_file_id: str, study_id: str, **kwargs):
        GeneExpression.__init__(self, mtx_path, study_file_id, study_id)

        genes_path = kwargs.pop("gene_file")
        genes_ingest_file = IngestFiles(genes_path, self.ALLOWED_FILE_TYPES)
        self.genes_file, genes_local_path = genes_ingest_file.resolve_path(
            genes_path)

        barcodes_path = kwargs.pop("barcode_file")
        barcodes_ingest_file = IngestFiles(
            barcodes_path, self.ALLOWED_FILE_TYPES)
        self.barcodes_file, barcodes_local_path = barcodes_ingest_file.resolve_path(
            barcodes_path
        )

        mtx_ingest_file = IngestFiles(mtx_path, self.ALLOWED_FILE_TYPES)
        self.mtx_file, self.mtx_local_path = mtx_ingest_file.resolve_path(
            mtx_path)
        # A list ['N', 'K', 'M'] that represents a gene-barcode matrix where N
        # is the gene index, M is the barcode index, and K is the expresion score
        # for the given gene index
        self.mtx_description = linecache.getline(
            self.mtx_local_path, 2).split()

    def execute_ingest(self):
        # import pdb
        # pdb.set_trace()
        self.extract_feature_barcode_matrices()
        # import pdb
        # pdb.set_trace()
        yield from self.transform()
        self.close()

    def extract_mtx(self, value):
        """
        Zgreps by gene index from mtx file to enhance performance/scale
        """
        return subprocess.run(['zgrep', f'^{value}\s', self.mtx_local_path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE
                              ).stdout.decode('utf-8')

    def extract_feature_barcode_matrices(self):
        """
        Sets relevant iterables for the gene and barcode file of the MTX bundle
        """
        self.genes = [g.strip().strip('\"')
                      for g in self.genes_file.readlines()]
        self.cells = [c.strip().strip('\"')
                      for c in self.barcodes_file.readlines()]

    def transform(self):
        """
        Transforms bundle into data models.
            - yields 5 gene models and its corresponding data array models
                at a time.
            - Tracks amount of time it takes to create models.
        Returns:
            tuple (generator):
                gene_models (list): gene model documents
                data_arrays (list): data arrays that correspond to gene models
        """
        num_processed = 0
        start_time = datetime.datetime.now()
        gene_models = []
        data_arrays = []
        # Create gene model for all cells available in file
        for all_cell_model in self.set_data_array_cells(
                self.cells, ObjectId()):
            data_arrays.append(all_cell_model)
        for mtx_gene_idx in range(int(self.mtx_description[0]) - 1):
            exp_scores = []
            exp_cells = []
            id = ObjectId()
            # Extract rows:str from mtx file that match current mtx_gene_idx
            # matched_rows = [[N, K, M],[N, K, M],...N] where N = mtx_gene_idx
            # Rows in mtx file are 1 based and represented as strings => have
            # to add 1 to mtx_gene_idx and convert to string.
            matched_rows = self.extract_mtx(
                str(mtx_gene_idx + 1)).split("\n")[:-1]
            # Grab gene_id and gene from features file
            # Location of gene in features file = mtx_gene_idx + 1
            gene_id, gene = self.genes[int(mtx_gene_idx)].split('\t')
            gene_models.append(self.Model(
                {
                    'name': gene,
                    'searchable_name': gene.lower(),
                    'study_file_id': self.study_file_id,
                    'study_id': self.study_id,
                    'gene_id': gene_id,
                    '_id': id,
                }
            ))
            # Collect all cells and expression scores associated with a gene
            for row in matched_rows:
                raw_gene_idx, raw_barcode_idx, raw_exp_score = row.split()
                exp_cell = self.cells[int(raw_barcode_idx) - 1]
                exp_score = round(float(raw_exp_score), 3)
                exp_cells.append(exp_cell)
                exp_scores.append(exp_score)
            for cell_data_array in self.set_data_array_gene_cell_names(
                    gene,
                    id,
                    exp_cells):
                data_arrays.append(cell_data_array)
            for gene_data_array in self.set_data_array_gene_expression_values(gene, id, exp_scores):
                data_arrays.append(gene_data_array)
            if len(gene_models) > 5:
                num_processed += len(gene_models)
                yield (gene_models, data_arrays)
                print(f'Processed {num_processed} models, {str(datetime.datetime.now() - start_time)} elapsed')
                gene_models = []
                data_arrays = []
        yield (gene_models, data_arrays)
        num_processed += len(gene_models)
        print(f'Processed {num_processed} models, {str(datetime.datetime.now() - start_time)}')
        gene_models = []
        data_arrays = []

    @trace
    def close(self):
        """
        Closes file barcode and features file
        """
        self.genes_file.close()
        self.barcodes_file.close()
