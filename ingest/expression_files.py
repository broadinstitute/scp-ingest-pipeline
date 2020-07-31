"""A base class for expression files.

DESCRIPTION
This base class provides functions to create dataArrays for gene expression
files.

PREREQUISITES
Must have python 3.6 or higher.
"""
import abc
from dataclasses import dataclass
from mypy_extensions import TypedDict
from typing import List  # noqa: F401
import ntpath
from bson.objectid import ObjectId

try:
    from ingest_files import DataArray
    from monitor import setup_logger
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import DataArray
    from .monitor import setup_logger


class GeneExpression:
    __metaclass__ = abc.ABCMeta
    COLLECTION_NAME = "genes"
    info_logger = setup_logger(__name__, "info.txt")

    @dataclass
    class Model(TypedDict):
        name: str
        # downcase version of 'name'
        searchable_name: str
        study_file_id: ObjectId
        study_id: ObjectId
        _id: ObjectId
        gene_id: str = None

    def __init__(self, file_path: str, study_id: str, study_file_id: str):
        self.study_id = ObjectId(study_id)
        self.study_file_id = ObjectId(study_file_id)
        self.head, self.tail = ntpath.split(file_path)
        self.cluster_name = self.tail or ntpath.basename(self.head)
        self.extra_log_params = {"study_id": self.study_id, "duration": None}

    def create_gene_model(self, gene_name: str, study_file_id, study_id,
                          gene_id: int, model_id: int):
        return self.Model(
            {
                'name': gene_name,
                'searchable_name': gene_name.lower(),
                'study_file_id': study_file_id,
                'study_id': study_id,
                'gene_id': gene_id,
                '_id': model_id,
            }
        )

    @abc.abstractmethod
    def transform(self):
        """Abstract method for transforming expression data into Gene data model"""

    def set_data_array(self, name, array_type, values: List, linear_data_type,linear_data_id):
        """
        Sets annotation array
        """

        for model in DataArray(
            name,
            self.cluster_name,
            array_type,
            values,
            linear_data_type,
            linear_data_id,
            self.study_id,
            self.study_file_id,
        ).get_data_array():
            yield model

    def set_data_array_gene_cell_names(
        self, name: str, linear_data_id: str, values: List
    ):
        """Sets DataArray for cell names associated to a single gene. This
        DataArray contains cell names that had significant (i.e. non-zero)
        expression for a gene.
        """
        for model in DataArray(
            f"{name} Cells",
            self.cluster_name,
            "cells",
            values,
            "Gene",
            linear_data_id,
            self.study_id,
            self.study_file_id,
        ).get_data_array():
            yield model

    def set_data_array_gene_expression_values(
        self, name: str, linear_data_id: str, values: List
    ):
        """
        Sets DataArray for expression values for a gene. This is an array of
        significant (i.e. non-zero) expression values for a gene.
        """
        for model in DataArray(
            f"{name} Expression",
            self.cluster_name,
            "expression",
            values,
            "Gene",
            linear_data_id,
            self.study_id,
            self.study_file_id,
        ).get_data_array():
            yield model
