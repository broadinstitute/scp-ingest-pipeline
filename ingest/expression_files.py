"""
An abstract Base class for expression files.

DESCRIPTION
This base class provides functions to create dataArrays for gene expression
files.

PREREQUISITES
Must have python 3.6 or higher.
"""
from ingest_files import IngestFiles, DataArray
import abc
from dataclasses import dataclass
from mypy_extensions import TypedDict
from typing import List  # noqa: F401
import ntpath


class GeneExpression(IngestFiles):
    __metaclass__ = abc.ABCMeta
    COLLECTION_NAME = 'genes'

    @dataclass
    class Model(TypedDict):
        name: str
        # downcase version of 'name'
        searchable_name: str
        study_file_id: str
        study_id: str
        gene_id: str = None

    def __init__(
        self,
        file_path: str,
        study_id: str,
        study_file_id: str,
        allowed_file_types: str = None,
        open_as=None,
    ):
        self.study_id = study_id
        self.study_file_id = study_file_id
        self.head, self.tail = ntpath.split(file_path)
        self.cluster_name = self.tail or ntpath.basename(self.head)
        if open_as is not None:
            IngestFiles.__init__(
                self, file_path, allowed_file_types, open_as='dataframe'
            )

    @abc.abstractmethod
    def transform(self):
        """Abstract method for transforming expression data into Gene data model"""

    # This will end up being a class method
    @abc.abstractmethod
    def set_dataArray(self):
        """An abstract method that will be implemented by inherrited classes.
        Each expression file will have its own implementation of setting the
        DataArray.
        """

    def set_data_array_cells(self, values, linear_data_id):
        """Sets DataArray for cells that were observed in a
        expression matrix."""
        return DataArray(
            {
                'name': f'{self.cluster_name} Cells',
                'cluster_name': self.cluster_name,
                'array_type': 'cells',
                'linear_data_type': 'Study',
                'linear_data_id': linear_data_id,
                'vaules': values,
            }
        )

    def set_data_array_gene_cell_names(self, name, linear_data_id, values):
        """Sets DataArray for gene cell names. This DataArray contains cell
        names that had significant (i.e. non-zero) expression for a gene. """
        return DataArray(
            {
                'name': f'{name} Cells',
                'cluster_name': self.cluster_name,
                'array_type': 'cells',
                'linear_data_type': 'Gene',
                'linear_data_id': linear_data_id,
                'values': values,
            }
        )

    def set_data_array_gene_expression_values(self, name, linear_data_id, values):
        """ Sets DataArray for gene expression values. This is an array of
        significant (i.e. non-zero) expression values for a gene. """
        return DataArray(
            {
                'name': f'{name} Expression',
                'cluster_name': self.cluster_name,
                'array_type': 'expression',
                'linear_data_type': 'for Gene',
                'linear_data_id': linear_data_id,
                'value': values,
            }
        )
