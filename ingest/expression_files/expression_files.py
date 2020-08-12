"""A base class for expression files.

DESCRIPTION
This base class provides functions to create dataArrays for gene expression
files.

PREREQUISITES
Must have python 3.6 or higher.
"""
import abc
import datetime
import ntpath
import sys
from dataclasses import dataclass
from typing import List, Dict  # noqa: F401

from bson.objectid import ObjectId
from mypy_extensions import TypedDict
from pymongo.errors import BulkWriteError

try:
    sys.path.append("..")
    # Used when importing as external package, e.g. imports in single_cell_portal cod
    from ingest_files import DataArray
    from monitor import setup_logger
    from mongo_connection import MongoConnection
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..ingest_files import DataArray
    from ..monitor import setup_logger
    from ..mongo_connection import MongoConnection


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
        head, tail = ntpath.split(file_path)
        self.cluster_name = tail or ntpath.basename(head)
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

    @staticmethod
    def check_unique_cells(cell_names: List, study_id, client):
        """Checks cell names against database to confirm matrix contains unique
            cell names

         Parameters:
            cell_names (List[str]): List of cell names in matrix
            study_id (ObjectId): The study id the cell names belong to
            client: MongoDB client
        """
        COLLECTION_NAME = "data_arrays"
        query = {
            "$and": [
                {"linear_data_type": "Study"},
                {"array_type": "cells"},
                {"study_id": study_id},
            ]
        }
        # Returned fields from query results
        field_names = {"values": 1, "_id": 0}
        # Dict = {values_1: [<cell names>]... values_n:[<cell names>]}
        query_results: List[Dict] = list(
            client[COLLECTION_NAME].find(query, field_names)
        )
        # Query did not return results
        if not query_results:
            return True
        # Flatten query results
        existing_cells = [
            values
            for cell_values in query_results
            for values in cell_values.get("values")
        ]
        dupes = set(existing_cells) & set(cell_names)
        if len(dupes) > 0:
            error_string = f'Expression file contains {len(dupes)} cells that also exist in another expression file. '
            # add the first 3 duplicates to the error message
            error_string += f'Duplicates include {", ".join(list(dupes)[:3])}'
            raise ValueError(error_string)
        return True


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

    def load(self, gene_docs: List, data_array_docs: List):
        start_time = datetime.datetime.now()
        self.insert(gene_docs, self.COLLECTION_NAME)
        self.insert(data_array_docs, "data_arrays")
        self.info_logger.info(
            f"Time to load {len(gene_docs) + len(data_array_docs)} models: {str(datetime.datetime.now() - start_time)}"
        )

    def insert(self, docs: List, collection_name: str):
        try:
            self.mongo_connection._client[collection_name].insert_many(
                docs, ordered=False
            )
        except BulkWriteError as bwe:
            raise BulkWriteError(
                f"Error caused by inserting into collection '{collection_name}': {bwe.details}"
            )
        except Exception as e:
            raise Exception(
                f"Error caused by inserting into collection '{collection_name}': {e}"
            )
