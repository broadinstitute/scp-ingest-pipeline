"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""

import collections
from typing import List  # noqa: F401

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


class Dense(GeneExpression, IngestFiles):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        self.tracer = kwargs.pop("tracer")
        GeneExpression.__init__(self, file_path, study_file_id, study_id)
        IngestFiles.__init__(
            self, file_path, allowed_file_types=self.ALLOWED_FILE_TYPES
        )
        self.matrix_params = kwargs
        self.gene_names = []
        csv_file, self.open_file_object = self.open_file(self.file_path)
        self.header = next(csv_file)

    def validate_unique_header(self):
        """Validates header has no duplicate values"""
        if len(set(self.header)) != len(self.header):
            self.error_logger.error(
                "Duplicated header values are not allowed", extra=self.extra_log_params
            )
            return False
        return True

    def validate_gene_keyword(self):
        """Validates that 'Gene' is the first value in header"""
        # File is an R formatted file
        if (self.header[-1] == '') and (self.header[0].upper() != 'GENE'):
            pass
        else:
            # If not R formatted file then first cell must be 'GENE'
            if self.header[0].upper() != 'GENE':
                self.error_logger.error(
                    f'First cell in header must be GENE not {self.header[0]}',
                    extra=self.extra_log_params,
                )
                return False
        return True

    @trace
    def preprocess(self):
        """Determines if file is R-formatted. Creates dataframe (df)"""
        csv_file, open_file_object = self.open_file(self.file_path)
        dtypes = {'GENE': object}

        # # Remove white spaces and quotes
        header = [col_name.strip().strip('\"') for col_name in self.header]
        # See if R formatted file
        if (header[-1] == '') and (header[0].upper() != 'GENE'):
            header.insert(0, 'GENE')
            # Although the last column in the header is empty, python treats it
            # as an empty string,[ ..., ""]
            header = header[0:-1]
        else:
            header[0] = header[0].upper()
        # Set dtype for expression values to floats
        dtypes.update({cell_name: 'float' for cell_name in header[1:]})
        self.df = self.open_file(
            self.file_path,
            open_as='dataframe',
            names=header,
            skiprows=1,
            dtype=dtypes,
            # chunksize=100000, Save for when we chunk data
        )[0]

    @trace
    def transform(self):
        """Transforms dense matrix into gene data model.
        """
        # Holds gene name and gene model for a single gene
        GeneModel = collections.namedtuple('Gene', ['gene_name', 'gene_model'])
        for gene in self.df['GENE']:
            formatted_gene_name = gene.strip().strip('\"')
            self.info_logger.info(
                f'Transforming gene :{gene}', extra=self.extra_log_params
            )
            if gene in self.gene_names:
                raise ValueError(f'Duplicate gene: {gene}')
            self.gene_names.append(gene)
            id = ObjectId()

            yield GeneModel(
                # Name of gene as observed in file
                gene,
                self.Model(
                    {
                        'name': formatted_gene_name,
                        'searchable_name': formatted_gene_name.lower(),
                        'study_file_id': self.study_file_id,
                        'study_id': self.study_id,
                        '_id': id,
                        'gene_id': self.matrix_params['gene_id']
                        if 'gene_id' in self.matrix_params
                        else None,
                    }
                ),
            )

    @trace
    def set_data_array(
        self,
        linear_data_id,
        unformatted_gene_name,
        gene_name,
        create_cell_data_array=False,
    ):
        """Creates data array for expression, gene, and cell values."""
        input_args = locals()
        gene_df = self.df[self.df['GENE'] == unformatted_gene_name]
        # Get list of cell names
        cells = self.df.columns.tolist()[1:]
        if create_cell_data_array:
            self.info_logger.info(
                f'Creating cell data array for gene : {gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_cells(cells, linear_data_id)

        # Get row of expression values for gene
        # Round expression values to 3 decimal points
        cells_and_expression_vals = gene_df[cells].round(3).to_dict('records')[0]
        # Filter out expression values = 0
        cells_and_expression_vals = dict(
            filter(lambda k_v: k_v[1] > 0, cells_and_expression_vals.items())
        )
        observed_cells = list(cells_and_expression_vals)
        observed_values = list(cells_and_expression_vals.values())

        # Return significant data (skip if no expression was observed)
        if len(observed_values) > 0:
            self.info_logger.info(
                f'Creating cell names data array for gene: {unformatted_gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_gene_cell_names(
                unformatted_gene_name, linear_data_id, observed_cells
            )
            self.info_logger.info(
                f'Creating gene expression data array for gene: {unformatted_gene_name}',
                extra=self.extra_log_params,
            )
            yield from self.set_data_array_gene_expression_values(
                unformatted_gene_name, linear_data_id, observed_values
            )

    def validate_format(self):
        return all([self.validate_unique_header(), self.validate_gene_keyword()])

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.open_file_object.close()
