"""A base class for expression files.

DESCRIPTION
This base class provides functions to create dataArrays for gene expression
files.

PREREQUISITES
Must have python 3.6 or higher.
"""
import abc
import ntpath
import types
from dataclasses import dataclass
from typing import List  # noqa: F401

from bson.objectid import ObjectId
from mypy_extensions import TypedDict
import datetime

from pymongo import InsertOne, MongoClient
from pymongo.errors import BulkWriteError

try:
    from ingest_files import DataArray
    from monitor import setup_logger
    from connection import MongoConnection
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import DataArray
    from .monitor import setup_logger
    from connection import MongoConnection


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
        self.mongo_connection = MongoConnection()

    @abc.abstractmethod
    def transform(self):
        """Abstract method for transforming expression data into Gene data model"""

    @abc.abstractmethod
    def set_dataArray(self):
        """An abstract method that will be implemented by inherrited classes.
        Each expression file will have its own implementation of setting the
        DataArray with expression data."""

    def set_data_array_cells(self, values: List, linear_data_id):
        """Sets DataArray for cells that were observed in an
        expression matrix."""
        for model in DataArray(
            f"{self.cluster_name} Cells",
            self.cluster_name,
            "cells",
            values,
            "Study",
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
        expression for a gene. """

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
        """ Sets DataArray for expression values for a gene. This is an array of
        significant (i.e. non-zero) expression values for a gene. """

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

    def load_expression_file(self, gene_docs: List, data_array_documents: List):
        """
        
        """
        gene_doc_bulk_write_results = None
        data_array_bulk_write_results = None
        start_time = datetime.datetime.now()

        # Creating Mongo bulk operations
        data_array_bulk_operations = list(
            map(lambda model: InsertOne(model), data_array_documents))
        gene_model_bulk_operations = list(
            map(lambda model: InsertOne(model), gene_docs))

        # Try writing data_array_colection
        try:
            self.mongo_connection.client['data_arrays'].bulk_write(
                data_array_bulk_operations,  ordered=False
            )
        except BulkWriteError as bwe:
            print(f"error caused by data docs : {bwe.details}")
            self.error_logger.error(bwe.details, extra=self.extra_log_params)
            return False

        except Exception as e:
            print(f"error caused by data docs : {e}")
            self.error_logger.error(e, extra=self.extra_log_params)
            return False

        # Try writing gene docs
        try:
            gene_doc_bulk_write_results = self.mongo_connection.client[self.COLLECTION_NAME].bulk_write(
                gene_model_bulk_operations,  ordered=False
            )
        except BulkWriteError as bwe:
            print(f"error caused by gene docs : {bwe.details}")
            self.error_logger.error(bwe.details, extra=self.extra_log_params)
            return False

        except Exception as e:
            print(f"error caused by gene docs : {e}")
            self.error_logger.error(e, extra=self.extra_log_params)
            return False
        print(f'Time to load {len(data_array_bulk_operations) + len(gene_model_bulk_operations)} models: {str(datetime.datetime.now() - start_time)}')
        return True
