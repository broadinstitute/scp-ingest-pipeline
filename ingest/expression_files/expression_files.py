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
import copy
from dataclasses import dataclass
from typing import List, Dict, Generator  # noqa: F401

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
    DATA_ARRAY_BATCH_SIZE = 1000
    # Logger provides more details
    dev_logger = setup_logger(__name__, "log.txt", format="support_configs")

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
        self.mongo_connection = MongoConnection()
        # Common data array kwargs
        self.data_array_kwargs = {
            "cluster_name": self.cluster_name,
            "study_file_id": self.study_file_id,
            "study_id": self.study_id,
        }

    @abc.abstractmethod
    def transform(self):
        """Abstract method for transforming expression data into data models."""

    @abc.abstractmethod
    def execute_ingest(self):
        """Abstract method for parsing expression data into MongoDB."""

    @staticmethod
    @abc.abstractmethod
    def check_valid():
        """Abstract method for validating expression matrices."""

    @staticmethod
    def create_gene_model(
        *ignore, name: str, study_file_id, study_id, _id: int, gene_id: str = None
    ):
        """Creates a gene model for a single gene.
            This function accepts keyword arguments only. An error will be raise
                when positional or additional keyword arguments are passed in.
        """
        if ignore:
            raise TypeError("Position arguments are not accepted.")
        return GeneExpression.Model(
            {
                "name": name,
                "searchable_name": name.lower(),
                "study_file_id": study_file_id,
                "study_id": study_id,
                "gene_id": gene_id,
                "_id": _id,
            }
        )

    @staticmethod
    def check_unique_cells(cell_names: List, study_id, client):
        """Checks cell names against database to confirm matrix contains unique
            cell names.

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
            ],
            "$nor": [{"name": "All Cells"}],
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
            error_string = (
                f"Expression file contains {len(dupes)} cells "
                "that also exist in another expression file."
            )

            # add the first 3 duplicates to the error message
            error_string += f'Duplicates include {", ".join(list(dupes)[:3])}'
            raise ValueError(error_string)
        return True

    @staticmethod
    def create_data_arrays(
        *ignore,
        # keyword arguments
        name: str,
        cluster_name: str,
        array_type: str,
        values: List,
        linear_data_type: str,
        linear_data_id,
        study_id,
        study_file_id,
    ) -> Generator:
        """
        Sets data array for expression data.
        This function accepts keyword arguments only. An error will be raise
                when positional or additional keyword arguments are passed in.
        """
        fn_kwargs = copy.copy(locals())
        # Positional arguments passed in
        if ignore:
            raise TypeError("Positional arguments are not accepted.")
        del fn_kwargs["ignore"]
        for model in DataArray(**fn_kwargs).get_data_arrays():
            yield model

    @staticmethod
    def insert(docs: List, collection_name: str, client):
        try:
            client[collection_name].insert_many(docs, ordered=False)
        except BulkWriteError as bwe:
            raise BulkWriteError(
                f"Error caused by inserting into collection '{collection_name}': {bwe.details}"
            )
        except Exception as e:
            raise Exception(
                f"Error caused by inserting into collection '{collection_name}': {e}"
            )

    def load(self, gene_docs: List, data_array_docs: List):
        start_time = datetime.datetime.now()
        GeneExpression.insert(gene_docs, self.COLLECTION_NAME, self.mongo_connection)
        GeneExpression.insert(data_array_docs, "data_arrays", self.mongo_connection)
        GeneExpression.dev_logger.info(
            f"Time to load {len(gene_docs) + len(data_array_docs)} models: {str(datetime.datetime.now() - start_time)}"
        )
