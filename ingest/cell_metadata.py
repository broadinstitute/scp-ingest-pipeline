"""Module for ingesting cell metadata files
DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.
PREREQUISITES
Must have python 3.6 or higher.
"""
import copy
from collections import defaultdict
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

from ingest_files import IngestFiles


class CellMetadata(IngestFiles):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]

    def __init__(self, file_path, file_id: str, study_accession: str, *args, **kwargs):

        IngestFiles.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, open_as="dataframe"
        )
        self.headers = self.get_next_line(increase_line_count=False)
        self.metadata_types = self.get_next_line(increase_line_count=False)
        # unique values for group-based annotations
        self.unique_values = []
        self.cell_names = []
        self.annotation_type = ["group", "numeric"]
        self.top_level_doc = self.create_documents(file_id, study_accession)
        self.data_subcollection = self.create_subdocuments()
        self.errors = defaultdict(list)
        self.ontology = defaultdict(lambda: defaultdict(set))
        self.type = defaultdict(list)
        self.cells = []

    def transform(self, row: List[str]) -> None:
        """ Add data from cell metadata files into data model"""
        for idx, column in enumerate(row):
            if idx != 0:
                # if annotation is numeric convert from string to float
                if self.metadata_types[idx].lower() == "numeric":
                    column = round(float(column), 3)
                elif self.metadata_types[idx].lower() == "group":
                    # Check for unique values
                    if column not in self.unique_values:
                        self.unique_values.append(column)
                # Get annotation name from header
                annotation = self.headers[idx]
                self.data_subcollection[annotation]["values"].append(column)
            else:
                # If column isn't an annotation value, it's a cell name
                self.cell_names.append(column)

    def create_documents(self, file_id, study_accession):
        """Creates top level documents for Cell Metadata data structure"""
        documents = {}

        # Each annotation value has a top level document
        for idx, value in enumerate(self.headers[1:]):
            # Copy document model so memory references are different
            copy_of_doc_model = copy.copy(
                {
                    "name": value,
                    "study_accession": study_accession,
                    "unique_values": [],
                    "annotation_type": self.metadata_types[idx + 1],
                    "file_id": file_id,
                }
            )
            documents[value] = copy_of_doc_model
        return documents

    def create_subdocuments(self):
        """Creates subdocuments for each annotation """
        sub_documents = {}
        for value in self.headers[1:]:
            # Copy subdocument model so memory references are different
            copy_of_subdoc_model = copy.copy(
                {"cell_names": self.cell_names, "values": []}
            )
            sub_documents[value] = copy_of_subdoc_model
        return sub_documents

    def get_collection_name(self):
        """Returns collection name"""
        return "cell_metadata"

    def get_subcollection_name(self):
        """Returns sub-collection name"""
        return "data"

    def validate_header_keyword(self):
        """Check metadata header row starts with NAME (case-insensitive).

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.headers[0].casefold() == "NAME".casefold():
            valid = True
            if self.headers[0] != "NAME":
                # ToDO - capture warning below in error report
                print(
                    f'Warning: metadata file keyword "NAME" provided as '
                    "{self.headers[0]}"
                )
        else:
            # line below and similar in next method have autoformat oddities
            self.errors["format"].append(
                "Error: Metadata file header row malformed, missing NAME"
            )
        return valid

    def validate_unique_header(self):
        """Check all metadata header names are unique.

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if len(self.headers[1:]) == len(set(self.headers[1:])):
            valid = True
        else:
            self.errors["format"].append(
                "Error:  Duplicate column headers in metadata file"
            )
        return valid

    def validate_type_keyword(self):
        """Check metadata second row starts with TYPE (case-insensitive).

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.metadata_types[0].casefold() == "TYPE".casefold():
            valid = True
            if self.metadata_types[0] != "TYPE":
                # ToDO - capture warning below in error report
                # investigate f-string formatting here
                print(
                    'Warning: metadata file keyword "TYPE" provided as '
                    "{self.metadata_types[0]}"
                )
        else:
            # check black autoformatting on this long line
            self.errors["format"].append(
                "Error:  Metadata file TYPE row malformed, missing TYPE"
            )
        return valid

    def validate_type_annotations(self):
        """Check metadata second row contains only 'group' or 'numeric'.

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        annot_err = False
        annots = []
        # skipping the TYPE keyword, iterate through the types
        # collecting invalid type annotations in list annots
        for t in self.metadata_types[1:]:
            if t not in self.annotation_type:
                # if the value is a blank space, store a higher visibility
                # string for error reporting
                if not t:
                    annots.append("<empty value>")
                else:
                    annots.append(t)
                annot_err = True
        if annot_err:
            self.errors["format"].append(
                (
                    'Error: TYPE declarations should be "group" or "numeric"; '
                    f'Invalid type(s): {", ".join(map(str, annots))}'
                )
            )
        else:
            valid = True
        return valid

    def validate_against_header_count(self):
        """Metadata header and type counts should match.

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if not len(self.headers) == len(self.metadata_types):
            self.errors["format"].append(
                "Error: {len(self.metadata_types)} TYPE declarations "
                f"for {len(self.headers)} column headers"
            )
        else:
            valid = True
        return valid

    def validate_format(self):
        """Check all metadata file format criteria for file validity
        """
        self.validate_header_keyword()
        self.validate_type_keyword()
        self.validate_type_annotations()
        self.validate_unique_header()
        self.validate_against_header_count()
        if self.errors["format"]:
            valid = False
        else:
            valid = True
        return valid
