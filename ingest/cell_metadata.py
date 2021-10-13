"""Module for ingesting cell metadata files

DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.

PREREQUISITES
Must have python 3.6 or higher.
"""
import collections
import ntpath
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
import copy
import json

from bson.objectid import ObjectId
from mypy_extensions import TypedDict

try:
    # Used when importing internally and in tests
    from annotations import Annotations
    from ingest_files import DataArray, IngestFiles
    from validation.validate_metadata import (
        report_issues,
        validate_input_metadata,
        write_metadata_to_bq,
    )
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .annotations import Annotations
    from .ingest_files import DataArray, IngestFiles


class CellMetadata(Annotations):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    COLLECTION_NAME = "cell_metadata"
    # File location for metadata json convention
    JSON_CONVENTION = (
        "../schema/alexandria_convention/alexandria_convention_schema.json"
    )

    def __init__(
        self,
        file_path: str,
        study_id: ObjectId,
        study_file_id: ObjectId,
        *args,
        **kwargs,
    ):

        self.study_accession = kwargs.pop("study_accession")
        Annotations.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, study_id, study_file_id
        )
        self.cell_names = []
        # lambda below initializes new key with nested dictionary as value and avoids KeyError
        self.issues = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.ontology = defaultdict(lambda: defaultdict(list))
        self.ontology_label = dict()
        self.cells = []
        self.numeric_array_columns = {}
        self.kwargs = kwargs

    # This model pertains to columns from cell metadata files
    @dataclass
    class Model(TypedDict):
        # value from column header
        name: str
        annotation_type: str
        study_file_id: ObjectId
        study_id: ObjectId
        # unique values from "group" type annotations
        values: List

    # More work needs to be done to fully remove ingest from IngestPipeline
    # Tracked in SCP-3023
    def execute_ingest(self):
        """ Method for ingesting cell metadata files."""
        for metadata_model in self.transform():
            yield metadata_model

    def update_numeric_array_columns(self, annotation_name):
        if not self.numeric_array_columns.values():
            self.numeric_array_columns[annotation_name] = True

    # Will evolve to do cross file validation
    def validate(self, validate_against_convention=False):
        """ Runs all validation checks """
        if validate_against_convention:
            return self.is_valid_format()
        else:
            return all([self.is_valid_format(), self.validate_numeric_annots()])

    def is_valid_format(self):
        """Validates format by calling all format validation methods"""
        return all(
            [self.validate_header_for_coordinate_values(), self.validate_format()]
        )

    def validate_header_for_coordinate_values(self):
        """Cell metadata files should not have coordinates in header
        :return: boolean True if coordinates are not in header, otherwise False
        """
        lower_cased_headers = [header.lower() for header in self.headers]
        valid = not any(
            [coordinate in ("x", "y", "z") for coordinate in lower_cased_headers]
        )
        if valid:
            return True
        else:
            msg = "Header names can not be coordinate values x, y, or z (case insensitive)"
            self.store_validation_issue("error", "format", msg)
            return False

    def conforms_to_metadata_convention(self):
        """ Determines if cell metadata file follows metadata convention"""
        convention_file_object = IngestFiles(self.JSON_CONVENTION, ["application/json"])
        json_file = convention_file_object.open_file(self.JSON_CONVENTION)
        convention = json.load(json_file)

        import_to_bq = self.kwargs["bq_dataset"] and self.kwargs["bq_table"]
        validate_input_metadata(self, convention, bq_json=import_to_bq)

        json_file.close()
        return not report_issues(self)

    def transform(self):
        """ Builds cell metadata model"""
        AnnotationModel = collections.namedtuple(
            "AnnotationModel", ["annot_header", "model"]
        )
        for annot_header in self.file.columns:
            annot_name = annot_header[0]
            annot_type = annot_header[1]
            # When file is conventional and contains numeric arrays
            # the annotation type is changed to group for visualization purposes
            stored_mongo_annot_type: str = (
                annot_type
                if not self.numeric_array_columns.get(annot_name)
                else "group"
            )
            yield AnnotationModel(
                annot_header,
                self.Model(
                    {
                        "name": annot_name,
                        "annotation_type": stored_mongo_annot_type,
                        # unique values from "group" type annotations else []
                        "values": list(self.file[annot_header].unique())
                        if annot_type == "group"
                        else [],
                        "study_file_id": self.study_file_id,
                        "study_id": self.study_id,
                    }
                ),
            )

    def set_data_array(self, linear_data_id: str, annot_header: str):

        data_array_attrs = copy.copy(locals())
        del data_array_attrs["annot_header"]
        del data_array_attrs["self"]
        annot_name = annot_header[0]
        head, tail = ntpath.split(self.file_path)
        base_data_array_model = {
            "cluster_name": tail or ntpath.basename(head),
            "values": list(self.file[annot_header]),
            "study_file_id": self.study_file_id,
            "study_id": self.study_id,
        }
        # This is an array (or group of arrays) of every cell
        if annot_name.lower() == "name":
            base_data_array_model.update(
                {
                    "name": "All Cells",
                    "array_type": "cells",
                    "linear_data_type": "Study",
                }
            )
        # data from cell metadata file that correspond to a column of data
        else:
            base_data_array_model.update(
                {
                    "name": annot_name,
                    "array_type": "annotations",
                    "linear_data_type": "CellMetadatum",
                }
            )
        return DataArray(**base_data_array_model, **data_array_attrs).get_data_arrays()

    def yield_by_row(self) -> None:
        """ Yield row from cell metadata file"""
        for row in self.file.itertuples(index=False):
            dict_row = row._asdict()
            yield dict_row.values()
