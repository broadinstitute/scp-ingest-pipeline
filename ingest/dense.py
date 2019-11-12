"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""

import collections
from typing import List  # noqa: F401

from expression_files import GeneExpression


class Dense(GeneExpression):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        GeneExpression.__init__(
            self,
            study_file_id,
            study_id,
            file_path=file_path,
            allowed_file_types=self.ALLOWED_FILE_TYPES,
            open_as='dataframe',
        )
        self.hasAll = True
        # Remove from dictionary any keys that have value=None
        self.matrix_params = kwargs
        # Remove trailing white spaces, and quotes from column names
        self.file.rename(
            columns=lambda x: x.strip().strip('\"').strip('\''), inplace=True
        )
        self.file.rename(
            columns=lambda x: x.strip().strip('\'').strip('\"'), inplace=True
        )

    def transform_expression_data_by_gene(self):
        """Transforms dense matrix into firestore data model for genes.

        Args:
            lines : List[str]
                Lines from dense matrix file

        Returns:
                transformed_data : List[Gene]
                A list of Gene objects
        """

        Gene_Model = collections.namedtuple('Gene', ['gene_name', 'gene_model'])
        for gene in self.file['GENE']:
            # Remove white spaces and quotes
            formatted_gene_name = gene.strip().strip('\"').strip('\'').strip('\"')
            yield Gene_Model(
                # Name of gene as observed in file
                gene,
                self.Model(
                    {
                        'name': formatted_gene_name,
                        'searchable_name': formatted_gene_name.lower(),
                        'study_file_id': self.study_file_id,
                        'study_id': self.study_id,
                        'gene_id': self.matrix_params['gene_id']
                        if 'gene_id' in self.matrix_params
                        else None,
                    }
                ),
            )

    def set_data_array(
        self,
        unformatted_gene_name,
        gene_name,
        linear_data_id,
        create_cell_DataArray=False,
    ):
        input_args = locals()
        gene_df = self.file[self.file['GENE'] == unformatted_gene_name]

        if create_cell_DataArray:
            yield self.set_data_array_cells(self.file['GENE'].tolist())
        else:
            # Get list of cell names
            cells = self.file.columns.tolist()[1:]
            # Get row of expression values for gene
            # Round expression values to 3 decimal points
            # Insure data type of float
            cells_and_expression_vals = (
                gene_df[cells].round(3).astype(float).to_dict('records')[0]
            )
            # Filter out expression values = 0
            cells_and_expression_vals = dict(
                filter(lambda k_v: k_v[1] > 0, cells_and_expression_vals.items())
            )
            yield self.set_data_array_gene_cell_names(
                gene_name, linear_data_id, list(cells_and_expression_vals.keys())
            )
            yield self.set_data_array_gene_expression_values(
                gene_name, linear_data_id, list(cells_and_expression_vals.values())
            )

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.file.close()
