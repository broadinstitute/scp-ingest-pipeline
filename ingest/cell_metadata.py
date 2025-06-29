"""Module for ingesting cell metadata files

DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.

PREREQUISITES
Must have python 3.6 or higher.
"""

import collections
import ntpath
from collections import defaultdict, OrderedDict
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
        self.ontology = defaultdict(lambda: defaultdict(list))
        self.ordered_ontology = defaultdict(list)
        self.ordered_labels = defaultdict(list)
        self.synonym_updates = defaultdict(lambda: defaultdict(str))
        self.cells = []
        self.numeric_array_columns = {}
        self.modalities = kwargs.get("has_modality", None)
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
        """Method for ingesting cell metadata files."""
        for metadata_model in self.transform():
            yield metadata_model

    def update_numeric_array_columns(self, annotation_name):
        if not self.numeric_array_columns.values():
            self.numeric_array_columns[annotation_name] = True

    # Will evolve to do cross file validation
    def validate(self, validate_against_convention=False):
        """Runs all validation checks"""
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
            msg = 'Header names can not be coordinate values "X", "Y", or "Z" (case insensitive).'
            self.store_validation_issue(
                "error", msg, "format:cap:metadata-no-coordinates"
            )
            return False

    @staticmethod
    def make_multiindex_name(modality):
        """From modality, generate column name in multi-index format"""
        has_modality_name = "has_" + modality
        multiindex_name = (has_modality_name, 'group')
        return multiindex_name

    @staticmethod
    def create_boolean_modality_metadatum(df, modality):
        """Translate presence of single modality to boolean for BigQuery"""
        # check for empty cells (aka. nan) or empty strings
        modality_multiindex = CellMetadata.make_multiindex_name(modality)
        no_modality_info = df[modality_multiindex].isna() | df[
            modality_multiindex
        ].str.len().eq(0)
        bool_name = modality + "_bool"
        bool_multiindex = CellMetadata.make_multiindex_name(bool_name)
        # store inverse of no_modality_info (ie. True = has modality info)
        df[bool_multiindex] = ~no_modality_info
        return df

    def hide_modality_metadatum(self):
        """Move non-boolean modality info out of self.file"""
        m_to_hide = []
        m_to_rename = {}
        for m in self.modalities:
            has_m = "has_" + m
            bool_name = has_m + "_bool"
            m_to_hide.append(CellMetadata.make_multiindex_name(m))
            m_to_rename[bool_name] = has_m
        self.modality_urls = self.file.filter(m_to_hide, axis=1)
        self.file.drop(m_to_hide, axis=1, inplace=True)
        self.file.rename(columns=m_to_rename, inplace=True)
        return

    def booleanize_modality_metadata(self):
        """Translate presence of modality data to boolean for BigQuery
        If no modality data, self.files is unchanged
        """
        if self.modalities is not None:
            df = copy.deepcopy(self.file)
            for m in self.modalities:
                df = CellMetadata.create_boolean_modality_metadatum(df, m)
            self.file = df
            self.hide_modality_metadatum()

    @staticmethod
    def restore_modality_metadata(cm):
        """Restore modality file data for MongoDB ingest"""
        if cm.modalities is not None:
            m_to_drop = []
            for m in cm.modalities:
                m_to_drop.append(CellMetadata.make_multiindex_name(m))
            cm.file.drop(m_to_drop, axis=1, inplace=True)
            cm.file = cm.file.join(cm.modality_urls)
        return cm.file

    def conforms_to_metadata_convention(self):
        """Determines if cell metadata file follows metadata convention"""
        convention_file_object = IngestFiles(self.JSON_CONVENTION, ["application/json"])
        json_file = convention_file_object.open_file(self.JSON_CONVENTION)
        convention = json.load(json_file)

        validate_input_metadata(self, convention)

        json_file.close()
        return not report_issues(self)

    def transform(self):
        """Builds cell metadata model"""
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

            group = (
                True
                if (annot_type == "group" or stored_mongo_annot_type == "group")
                else False
            )
            # should not store annotations with >200 unique values for viz
            # annot_header is the column of data, which includes name and type
            # large is any annotation with more than 200 + 2 unique values
            unique_values = list(self.file[annot_header].unique())
            large = True if len(unique_values) > 202 else False

            yield AnnotationModel(
                annot_header,
                self.Model(
                    {
                        "name": annot_name,
                        "annotation_type": stored_mongo_annot_type,
                        # unique values from "group" type annotations else []
                        "values": unique_values if (group and not large) else [],
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
        """Yield row from cell metadata file"""
        for row in self.file.itertuples(index=False):
            dict_row = row._asdict()
            yield dict_row.values()
