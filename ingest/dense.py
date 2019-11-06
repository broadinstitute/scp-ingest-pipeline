"""Module for ingesting dense matrix files

DESCRIPTION
Module provides extract and transforms function for gene expression data for
an dense matrix.

PREREQUISITES
Must have python 3.6 or higher.
"""
from typing import List  # noqa: F401
from ingest_files import IngestFiles
from dataclasses import dataclass

from mypy_extensions import TypedDict

import collections
import ntpath


class Dense(IngestFiles):
    LINEAR_DATA_TYPE = 'Gene'
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    # This model pertains to columns from cell metadata files
    @dataclass
    class Model(TypedDict):
        name: str
        # downcase version of 'name'
        searchable_name: str
        study_file_id: str
        study_id: str
        gene_id: str = None

    def __init__(self, file_path, study_file_id, study_id, **kwargs):
        IngestFiles.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, open_as='dataframe'
        )
        self.study_file_id = study_file_id
        self.study_id = study_id
        # Remove from dictionary any keys that have value=None
        self.matrix_params = kwargs
        # Remove trailing white spaces, and quotes from column names
        self.file.rename(
            columns=lambda x: x.strip().strip('\"').strip('\''), inplace=True
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
            formatted_gene = gene.strip().strip('\"').strip('\'')
            yield Gene_Model(
                gene,
                self.Model(
                    {
                        'name': formatted_gene,
                        'searchable_name': formatted_gene.lower(),
                        'study_file_id': self.study_file_id,
                        'study_id': self.study_id,
                        'gene_id': self.matrix_params['gene_id']
                        if 'gene_id' in self.matrix_params
                        else None,
                    }
                ),
            )

    def set_data_array(self, gene_name, name, linear_data_id):
        input_args = locals()
        cells = self.file.columns.tolist()[1:]
        values = self.file.loc[self.file['GENE'] == gene_name].values[1:]

        # values = [round(float(value), 3) if float(value)>0 for value in values]

        input_args['name'] = f'#{name} Cells'

        # return self.DataArray({**input_args, **base_data_array_model})

    def set_cell_data_array(self, name, linear_data_id):
        head, tail = ntpath.split(self.file_path)
        base_data_array_model = {
            'cluster_name': tail or ntpath.basename(head),
            'array_type': 'cells',
            'linear_data_type': 'Gene',
            'linear_data_id': linear_data_id,
        }
        return base_data_array_model

    def close(self):
        """Closes file

        Args:
            None

        Returns:
            None
        """
        self.file.close()
