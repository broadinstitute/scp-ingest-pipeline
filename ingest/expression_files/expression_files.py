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
import copy
from dataclasses import dataclass
from typing import List, Dict, Generator  # noqa: F401

from bson.objectid import ObjectId
from mypy_extensions import TypedDict

try:
    from ingest_files import DataArray
    from monitor import setup_logger
    from mongo_connection import MongoConnection, graceful_auto_reconnect
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..ingest_files import DataArray
    from ..monitor import setup_logger
    from ..mongo_connection import MongoConnection, graceful_auto_reconnect


class GeneExpression:
    __metaclass__ = abc.ABCMeta
    COLLECTION_NAME = "genes"
    DATA_ARRAY_BATCH_SIZE = 1_000
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
            This function accepts keyword arguments only. An error will be raised
                if positional or additional keyword arguments are passed in.
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
    def is_raw_count_file(study_id, study_file_id, client):
        "Checks if study file is a raw count matrix"
        COLLECTION_NAME = "study_files"
        QUERY = {"_id": study_file_id, "study_id": study_id}

        study_file_doc = list(client[COLLECTION_NAME].find(QUERY)).pop()
        # Name of embedded document that holds 'is_raw_count_files is named expression_file_info.
        # If study files does not have document expression_file_info
        # field, "is_raw_count_files", will not exist.:
        if "expression_file_info" in study_file_doc.keys():
            return study_file_doc["expression_file_info"]["is_raw_counts"]
        else:
            return False

    @staticmethod
    def query_cells(study_id, client, query_kwargs):
        QUERY = {
            "$and": [
                {"linear_data_type": "Study"},
                {"array_type": "cells"},
                {"study_id": study_id},
            ],
            "$nor": [{"name": "All Cells"}],
            **query_kwargs,
        }
        # Returned fields from query results
        FIELD_NAMES = {"values": 1, "_id": 0}
        COLLECTION_NAME = "data_arrays"

        return list(client[COLLECTION_NAME].find(QUERY, FIELD_NAMES))

    @staticmethod
    def get_cell_names_from_study_file_id(study_id, study_file_id, client):
        """Returns cell names of study files of the same type. However, the cell names
            from study file id that's passed into the function are not included."""

        additional_query_kwargs = {}
        query_results = []
        study_files_ids = GeneExpression.get_study_expression_file_ids(
            study_id, study_file_id, client
        )
        # There are no study file ids of the same type
        if not study_files_ids:
            return None

        # If there are study files of the same type create and add filters
        study_files_ids_list = [
            study_files_id["_id"] for study_files_id in study_files_ids
        ]
        additional_query_kwargs["study_file_id"] = {"$in": study_files_ids_list}

        # Dict = {values_1: [<cell names>]... values_n:[<cell names>]}
        query_results: List[Dict] = GeneExpression.query_cells(
            study_id, client, additional_query_kwargs
        )

        # there are no results so return an empty array
        if not query_results:
            existing_cells = []
        else:
            # Flatten query results
            existing_cells = [
                values
                for cell_values in query_results
                for values in cell_values.get("values")
            ]
        return existing_cells

    @staticmethod
    def check_unique_cells(cell_names: List, study_id, study_file_id, client):
        """Method checks for unique cells  by matrix type


         Parameters:
            cell_names (List[str]): List of cell names in matrix
            study_id (ObjectId): The study id the cell names belong to
            study_file_id (ObjectId): The file id the cell names belong to
            client: MongoDB client
        """
        # raise error if there are no cell names present
        if not cell_names:
            raise ValueError("There were no cell names found in the header row of this matrix")

        existing_cells = GeneExpression.get_cell_names_from_study_file_id(
            study_id, study_file_id, client
        )
        if existing_cells:
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
    def get_study_expression_file_ids(
        study_id, current_study_file_id, client
    ) -> List[Dict]:
        """Returns study file ids that are of the same type as 'current_study_file_id' in the study """

        COLLECTION_NAME = "study_files"
        field_names = {"_id": 1}
        QUERY = {
            "$and": [{"study_id": study_id}],
            "file_type": {"$in": ["Expression Matrix", "MM Coordinate Matrix"]},
            "$nor": [{"_id": current_study_file_id}],
        }
        is_raw_count_files = GeneExpression.is_raw_count_file(
            study_id, current_study_file_id, client
        )
        if is_raw_count_files:
            QUERY["$and"].append(
                {"expression_file_info.is_raw_counts": is_raw_count_files}
            )

        else:
            # Legacy data entries will not have the embedded document 'expression_file_info'
            QUERY["$and"].append(
                {
                    "$or": [
                        {"expression_file_info.is_raw_counts": is_raw_count_files},
                        {"expression_file_info": None},
                    ]
                }
            )
        query_results = list(client[COLLECTION_NAME].find(QUERY, field_names))
        return query_results

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
        This function accepts keyword arguments only. An error will be raised
                if positional or additional keyword arguments are passed in.
        """
        fn_kwargs = copy.copy(locals())
        # Positional arguments passed in
        if ignore:
            raise TypeError("Positional arguments are not accepted.")
        del fn_kwargs["ignore"]
        for model in DataArray(**fn_kwargs).get_data_arrays():
            yield model

    @staticmethod
    @graceful_auto_reconnect
    def insert(docs: List, collection_name: str, client):
        client[collection_name].insert_many(docs, ordered=False)

    def load(self, docs: List, collection_name: List):
        start_time = datetime.datetime.now()
        GeneExpression.insert(docs, collection_name, self.mongo_connection._client)
        GeneExpression.dev_logger.info(
            f"Time to load {len(docs)} models: {str(datetime.datetime.now() - start_time)}"
        )

    def create_models(
        self,
        exp_cells: List,
        exp_scores: List,
        gene: str,
        gene_id: str,
        gene_models: List,
        data_arrays: List,
        num_processed: int,
        force=False,
    ):
        """Creates models for a given gene and batches them for loading if
            necessary.

        After creating models, the amount of data arrays are checked to see if
        models should be loaded into the database. Data arrays and gene models
        will be loaded if the threshold, as specified in DATA_ARRAY_BATCH_SIZE,
        is met.

        Returns:
            data_arrays: List[str]
            gene_models: List[str]
            num_processed: int
                Amount of gene models created.
            """
        current_data_arrays = []
        start_time = datetime.datetime.now()
        model_id = ObjectId()

        if gene:
            gene_models.append(
                GeneExpression.create_gene_model(
                    name=gene,
                    study_file_id=self.study_file_id,
                    study_id=self.study_id,
                    gene_id=gene_id,
                    _id=model_id,
                )
            )
        # Make data array models for genes with expression data
        if len(exp_scores) > 0:
            # Data array model for cells
            for cell_data_array in GeneExpression.create_data_arrays(
                name=f"{gene} Cells",
                array_type="cells",
                values=exp_cells,
                linear_data_type="Gene",
                linear_data_id=model_id,
                **self.data_array_kwargs,
            ):
                current_data_arrays.append(cell_data_array)
            # Data array model for expression values
            for exp_value_data_array in GeneExpression.create_data_arrays(
                name=f"{gene} Expression",
                array_type="expression",
                values=exp_scores,
                linear_data_type="Gene",
                linear_data_id=model_id,
                **self.data_array_kwargs,
            ):
                current_data_arrays.append(exp_value_data_array)
        this_batch_size = len(data_arrays) + len(current_data_arrays)
        # Determine if models should be batched/loaded
        if this_batch_size >= GeneExpression.DATA_ARRAY_BATCH_SIZE or force:
            if force:
                # Add new data arrays
                data_arrays += current_data_arrays
                current_data_arrays.clear()
            if len(data_arrays) > 0:
                self.load(data_arrays, DataArray.COLLECTION_NAME)
            if len(gene_models) > 0:
                self.load(gene_models, GeneExpression.COLLECTION_NAME)
            num_processed += len(gene_models)
            GeneExpression.dev_logger.info(
                f"Processed {num_processed} genes. "
                f"{str(datetime.datetime.now() - start_time)} "
                f"elapsed"
            )
            gene_models.clear()
            data_arrays.clear()
        data_arrays += current_data_arrays

        return data_arrays, gene_models, num_processed
