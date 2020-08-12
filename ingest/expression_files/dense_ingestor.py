"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

"""
import datetime
import sys
import math
from typing import List  # noqa: F401

from bson.objectid import ObjectId

try:
    from expression_files import GeneExpression

    sys.path.append("..")
    from ingest_files import IngestFiles

except ImportError:
    # Used when importing as external package, e.g. imports in
    # single_cell_portal code
    from .expression_files import GeneExpression
    from ..ingest_files import IngestFiles


class DenseIngestor(GeneExpression, IngestFiles):
    # ToDo SCP-2635
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(self, file_path, study_file_id, study_id, **matrix_kwargs):
        GeneExpression.__init__(self, file_path, study_file_id, study_id)
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
        # To allow additional optional keyword arguments like gene_id
        self.matrix_params = matrix_kwargs
        self.csv_file_handler, self.file_handler = self.open_file(self.file_path)
        self.gene_names = {}
        self.header = DenseIngestor.process_header(next(self.csv_file_handler))

    def matches_file_type(file_type):
        return file_type == "dense"

    def execute_ingest(self):
        # Method can only be executed once due to
        # dependency on position in text file.
        # Row after header is needed for R format validation
        first_row = next(self.csv_file_handler)
        DenseIngestor.check_valid(
            self.header,
            first_row,
            query_params=(self.study_id, self.mongo_connection._client),
        )
        # Reset csv reader to first gene row
        self.csv_file_handler = self.open_file(self.file_path)[0]
        next(self.csv_file_handler)
        for gene_docs, data_array_documents in self.transform():
            self.load(gene_docs, data_array_documents)

    @staticmethod
    def check_valid(header, first_row, query_params):
      error_messages = []

      try:
          DenseIngestor.check_unique_header(header)
      except ValueError as v:
          error_messages.append(str(v))
      try:
          DenseIngestor.check_gene_keyword(header, first_row)
      except ValueError as v:
          error_messages.append(str(v))

      try:
          DenseIngestor.check_header_valid_values(header)
      except ValueError as v:
          error_messages.append(str(v))

      try:
          GeneExpression.check_unique_cells(header, *query_params)
      except ValueError as v:
          error_messages.append(str(v))

      if len(error_messages) > 0:
          raise ValueError('; '.join(error_messages))

      return True


    @staticmethod
    def format_gene_name(gene):
        return gene.strip().strip("'\",")

    @staticmethod
    def process_header(header):
        return [value.strip("'\",") for value in header]

    @staticmethod
    def process_row(row: str):
        """
        Performs pre-processing steps for a single row that converts gene
        scores into floats.
        """

        def convert_to_float(value: str):
            # Remove white spaces and quotes
            value = value.strip("'\",")
            # Convert string to float and round 3 places
            return round(float(value), 3)

        result = map(convert_to_float, row)
        return list(result)

    @staticmethod
    def filter_expression_scores(scores: List, cells: List):
        """
        Filters non-zero gene scores and corresponding cell names

        Returns:
            tuple (generator):
                valid_expression_scores (list): non-zero gene scores
                associated_cells (list): cells
        """
        associated_cells = []
        valid_expression_scores = []
        for idx, expression_score in enumerate(scores):
            try:
                if (
                    expression_score != "0"
                    and expression_score is not None
                    and not str(expression_score).isspace()
                    and expression_score != ""
                ):
                    # Can't evaluate strings for Nan values w/o breaking code
                    if expression_score != 0 and not math.isnan(
                        float(expression_score)
                    ):
                        valid_expression_scores.append(expression_score)
                        # add one to account for gene name in scores list
                        associated_cells.append(cells[idx + 1])
            except Exception:
                raise ValueError("Score '{expression_score}' is not valid")
        return valid_expression_scores, associated_cells

    @staticmethod
    def check_unique_header(header: List):
        """Confirms header has no duplicate values"""
        if len(set(header)) != len(header):
            raise ValueError("Duplicate header values are not allowed")
        return True

    @staticmethod
    def check_header_valid_values(header: List[str]):
        """Validates there are no empty header values"""
        for value in header:
            if value == "" or value.isspace():
                raise ValueError("Header values cannot be blank")
            if value.lower() == "nan":
                raise ValueError(f"{value} is not allowed as a header value")

        return True

    @staticmethod
    def check_gene_keyword(header: List, row: List):
        """Validates that 'Gene' is the first value in header

        Parameters:
            header (List[str]): Header of the dense matrix
            row (List): A single row from the dense matrix needed for R
            format validation
        """
        # Check if file is an R formatted file
        # An "R formatted" file has one less entry in the header
        # row than each successive row. Also, "GENE" will not appear in header
        if header[0].upper() != "GENE":
            length_of_next_line = len(row)
            if (length_of_next_line - 1) == len(header):
                return True
            else:
                raise ValueError("Required 'GENE' header is not present")
        else:
            return True

    def transform(self):
        """Transforms dense matrix into gene data model."""
        start_time = datetime.datetime.now()
        self.error_logger.info(
            "Starting run at " + str(start_time), extra=self.extra_log_params
        )
        num_processed = 0
        gene_models = []
        data_arrays = []
        for all_cell_model in self.set_data_array_cells(self.header[1:], ObjectId()):
            data_arrays.append(all_cell_model)
        # Represents row as a list
        for row in self.csv_file_handler:
            valid_expression_scores, cells = DenseIngestor.filter_expression_scores(
                row[1:], self.header
            )
            numeric_scores = DenseIngestor.process_row(valid_expression_scores)
            gene = row[0]
            if gene in self.gene_names:
                raise ValueError(f"Duplicate gene: {gene}")
            self.gene_names[gene] = True
            formatted_gene_name = DenseIngestor.format_gene_name(gene)
            id = ObjectId()
            gene_models.append(
                self.Model(
                    {
                        "name": formatted_gene_name,
                        "searchable_name": formatted_gene_name.lower(),
                        "study_file_id": self.study_file_id,
                        "study_id": self.study_id,
                        "_id": id,
                        "gene_id": None,
                    }
                )
            )
            if len(valid_expression_scores) > 0:
                for gene_cell_model in self.set_data_array_gene_cell_names(
                    gene, id, cells
                ):
                    data_arrays.append(gene_cell_model)
                for (
                    gene_expression_values
                ) in self.set_data_array_gene_expression_values(
                    gene, id, numeric_scores
                ):
                    data_arrays.append(gene_expression_values)
                if len(gene_models) == 5:
                    num_processed += len(gene_models)
                    self.info_logger.info(
                        f"Processed {num_processed} models, "
                        f"{str(datetime.datetime.now() - start_time)} elapsed",
                        extra=self.extra_log_params,
                    )
                    yield gene_models, data_arrays
                    gene_models = []
                    data_arrays = []
        yield gene_models, data_arrays
        num_processed += len(gene_models)
        self.info_logger.info(
            f"Processed {num_processed} models, {str(datetime.datetime.now() - start_time)} elapsed",
            extra=self.extra_log_params,
        )
