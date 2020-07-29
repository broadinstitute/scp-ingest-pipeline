"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

"""
import collections
import csv
import datetime
import sys
from typing import List  # noqa: F401

from bson.objectid import ObjectId

try:
    from .expression_files import GeneExpression
    sys.path.append("../ingest")
    from ingest_files import IngestFiles

except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .expression_files import GeneExpression
    sys.path.append("../ingest")
    from .ingest_files import IngestFiles


class DenseIngestor(GeneExpression, IngestFiles):
    ALLOWED_FILE_TYPES = ["text/csv",
                          "text/plain", "text/tab-separated-values"]

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        GeneExpression.__init__(self, file_path, study_file_id, study_id)
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
        self.matrix_params = kwargs
        self.csv_file_handler, self.file_handler = self.open_file(self.file_path)
        self.gene_names = {}
        self.header = next(self.csv_file_handler)

    def execute_ingest(self):
        # Row after header is needed for R format validation
        row = next(self.csv_file_handler)
        if DenseIngestor.is_valid_format(self.header, row):
            # Reset csv reader to first gene row
            self.csv_file_handler= self.open_file(self.file_path)[0]
            next(self.csv_file_handler)
            for gene_docs, data_array_documents in self.transform():
                self.load(gene_docs, data_array_documents)
        else:
            raise ValueError("Dense matrix has invalid format")

    def matches_file_type(file_type):
        return file_type == 'dense'

    @staticmethod
    def format_gene_name(gene):
        return gene.strip().strip('\"')

    @staticmethod
    def process_row(row: str):
        """
        Performs pre-processing steps a single row that converts gene scores
        into floats.
        """
        def convert_to_float(value: str):
            # Remove white spaces and quotes
            value = value.strip("'\",")
            # Convert string to float and round 3 places
            return round(float(value), 3)
        result = map(convert_to_float, row[1:])
        return list(result)

    @staticmethod
    def filter_expression_scores(scores: List, cells:List):
        """
        Filters non-zero gene scores and corresponding
        cell names

        Returns:
            tuple (generator):
                valid_expression_scores (list): non-zero gene scores
                associated_cells (list): cells
        """
        associated_cells = []
        valid_expression_scores = []
        for idx, expression_score in enumerate(scores, 1):
            if expression_score != 0:
                valid_expression_scores.append(expression_score)
                associated_cells.append(cells[idx])
        return valid_expression_scores, associated_cells

    @staticmethod
    def is_valid_format(header, row):
        return all(
            [
                DenseIngestor.has_unique_header(header),
                DenseIngestor.has_gene_keyword(header, row),
            ]
        )
    @staticmethod
    def has_unique_header(header):
        """Validates header has no duplicate values"""
        if len(set(header)) != len(header):
            # Logger will replace this
            print(
                "Duplicate header values are not allowed"
            )
            return False
        return True

    @staticmethod
    def has_gene_keyword(header: List, row:List):
        """
        Validates that 'Gene' is the first value in header

        Parameters:
            header (List[str]): Header of the dense matrix
            row (List): A single row from the dense matrix needed for R format validation
        """
        # Check if file is an R formatted file
        # An "R formatted" file has one less entry in the header
        # row than each successive row. Also, "GENE" will not appear in header
        if header[0].upper() != 'GENE':
            length_of_next_line = len(row)
            if (length_of_next_line-1) == len(header):
                return True
            else:
                return False
        else:
            return True

    def transform(self):
        """Transforms dense matrix into gene data model."""
        start_time = datetime.datetime.now()
        self.error_logger.info('Starting run at ' +
                               str(start_time), extra=self.extra_log_params)
        print('Starting run at ' + str(start_time))
        num_processed = 0
        gene_models = []
        data_arrays = []
        gene_id = None
        for all_cell_model in self.set_data_array_cells(
                self.header[1:], ObjectId()):
            data_arrays.append(all_cell_model)
        # Represents row as a list
        for row in self.csv_file_handler:
            numeric_scores = DenseIngestor.process_row(row)
            valid_expression_scores, cells = DenseIngestor.filter_expression_scores(
                numeric_scores, self.header)
            cells = []
            gene = row[0]
            if gene in self.gene_names:
                raise ValueError(f'Duplicate gene: {gene}')
            self.gene_names[gene] = True
            print(f'Transforming gene :{gene}')
            formatted_gene_name = DenseIngestor.format_gene_name(gene)
            id = ObjectId()
            if 'gene_id' in self.matrix_params:
                gene_id = self.matrix_params['gene_id']
            gene_models.append(self.Model(
                {
                    'name': formatted_gene_name,
                    'searchable_name': formatted_gene_name.lower(),
                    'study_file_id': self.study_file_id,
                    'study_id': self.study_id,
                    '_id': id,
                    'gene_id': gene_id
                }
            )
            )
            if len(valid_expression_scores) > 0:
                for gene_cell_model in self.set_data_array_gene_cell_names(gene, id, cells):
                    data_arrays.append(gene_cell_model)
                for gene_expression_values in self.set_data_array_gene_expression_values(gene, id, valid_expression_scores):
                    data_arrays.append(gene_expression_values)
                if len(gene_models) == 5:
                    num_processed += len(gene_models)
                    print(f'Processed {num_processed} models, {str(datetime.datetime.now() - start_time)} elapsed')
                    self.error_logger.info(f'Processed {num_processed} models, {str(datetime.datetime.now() - start_time)} elapsed', extra=self.extra_log_params)
                    yield (gene_models, data_arrays)
                    gene_models = []
                    data_arrays = []
        yield (gene_models, data_arrays)
        num_processed += len(gene_models)
        print(f'Processed {num_processed} models, {str(datetime.datetime.now() - start_time)}')
        gene_models = []
        data_arrays = []
