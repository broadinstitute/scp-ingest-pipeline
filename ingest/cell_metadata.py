"""Module for ingesting cell metadata files
DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.
PREREQUISITES
Must have python 3.6 or higher.
"""
from collections import defaultdict
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

from ingest_files import IngestFiles

DOCUMENT_LIMIT_BYTES = 1_048_576


class CellMetadata(IngestFiles):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    COLLECTION_NAME = "cell_metadata"
    SUBCOLLECTION_NAME = "data"

    def __init__(self, file_path, file_id: str, study_accession: str, *args, **kwargs):

        IngestFiles.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, open_as="dataframe"
        )
        self.is_valid_file = self.validate_format()
        self.preproccess()
        self.cell_names = []
        self.study_accession = study_accession
        self.file_id = file_id
        self.annotation_type = ['group', 'numeric']
        # lambda below initializes new key with nested dictionary as value and avoids KeyError
        self.issues = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.ontology = defaultdict(lambda: defaultdict(list))
        self.type = defaultdict(list)
        self.cells = []

    def preproccess(self):
        # Lowercase second level. Example: NUMeric -> numeric
        self.file.rename(
            columns=lambda col_name: col_name.lower(), level=1, inplace=True
        )
        self.file.xs("group", axis=1, level=1).astype(str)
        # Find numeric columns and round to 3 decimals places
        numeric_col_df = self.file.select_dtypes(include=["int64", "float64"]).columns
        self.file[numeric_col_df] = self.file[numeric_col_df].round(3).astype(float)

    def transform(self) -> None:
        """ Add data from cell metadata files into data model"""
        # first column is cell names, therefore skip
        for column in self.file.columns[1:]:
            col_name = column[0]
            column_type = column[1]
            yield (
                column_type,
                {
                    "name": col_name,
                    "study_accession": self.study_accession,
                    # save unique values for group type annotations
                    "unique_values": list(self.file[column].unique())
                    if column_type == "group"
                    else None,
                    "annotation_type": column_type,
                    "file_id": self.file_id,
                },
                {
                    "cell_names": list(self.file.iloc[:, 0]),
                    "values": list(self.file[column]),
                },
            )

    def chunk_subdocuments(self, doc_name, doc_path, subdoc, annot_type):
        """Partitions cell metadata subdocuments into storage sizes that are
            less than 1,048,576 bytes. Storage size calculation figures are derived from:
            # https://cloud.google.com/firestore/docs/storage-size

        Yields:
            Subdocuments that are under 1,048,576 bytes.
        """

        size_of_cell_names_field = 10 + 1  # "cell_names" is 10 characters
        size_of_values_field = 6 + 1  # "values" is 6 characters
        starting_sum = (
            +len(doc_name)
            + 1
            + len(doc_path)
            + 1
            + len(self.SUBCOLLECTION_NAME)
            + 1
            + len(self.COLLECTION_NAME)
            + 1
        )
        start_index = 0
        float_storage = 8
        sum = starting_sum

        # All cells names:[] that are in subdoc
        cell_names = subdoc["cell_names"]
        # All values:[] that are in subdoc
        values = subdoc["values"]

        for index, (cell_name, value) in enumerate(zip(cell_names, values)):

            cell_name_storage = len(cell_name) + 1 + size_of_cell_names_field

            # Check annotation type because float and string values have
            # different storage values
            if annot_type == "numeric":
                value_storage = size_of_values_field + float_storage
            else:
                value_storage = len(value) + 1 + size_of_values_field
            sum = sum + value_storage + cell_name_storage
            # Subtract 32 based off of firestore storage guidelines for strings
            # and documents
            # This and other storage size calculation figures are derived from:
            # https://cloud.google.com/firestore/docs/storage-size
            if (sum + 32) > DOCUMENT_LIMIT_BYTES or cell_name == cell_names[-1]:
                if cell_name == cell_names[-1]:
                    end_index = index
                else:
                    end_index = index - 1
                # TODO: This can turn into a logging statement
                # Please do not remove this. It's needed for testing
                print(f"{sum} , {index}, {start_index} , {end_index}")
                yield {
                    "cell_names": cell_names[start_index:end_index],
                    "values": values[start_index:end_index],
                }
                # Reset sum and add storage size at current index
                sum = starting_sum + cell_name_storage + value_storage
                start_index = index

    def store_validation_issue(self, type, category, msg, associated_info=None):
        """Store validation issues in proper arrangement
        :param type: type of issue (error or warn)
        :param category: issue category (format, jsonschema, ontology)
        :param msg: issue message
        :param value: list of IDs associated with the issue
        """
        if associated_info:
            self.issues[type][category][msg].extend(associated_info)
        else:
            self.issues[type][category][msg] = None

    def validate_header_keyword(self):
        """Check metadata header row starts with NAME (case-insensitive).
        :return: boolean   True if valid, False otherwise
        """

        """Check all metadata header names are unique.
        :return: boolean   True if valid, False otherwise
        """

        valid = False
        if self.file.columns[0][0].upper() == "NAME":
            valid = True
            if self.file.columns[0][0] != "NAME":
                # ToDO - capture warning below in error report
                msg = (
                    f'Warning: metadata file keyword "NAME" provided as '
                    f"{self.file.columns[0][0]}"
                )
                self.store_validation_issue('warn', 'format', msg)
        else:
            msg = 'Error: Metadata file header row malformed, missing NAME. (Case Sensitive)'
            self.store_validation_issue('error', 'format', msg)
        return valid

    def validate_unique_header(self):
        """Check all metadata header names are unique.
        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if len(set(self.file.columns.labels[0])) == len(self.file.columns.labels[0]):
            valid = True
        else:
            msg = 'Error: Duplicate column headers in metadata file'
            self.store_validation_issue('error', 'format', msg)
        return valid

    def validate_type_keyword(self):
        """Check metadata second row starts with TYPE (case-insensitive).
        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.file.columns[0][1].upper() == "TYPE":
            valid = True
            if self.file.columns[0][1] != "TYPE":
                # ToDO - capture warning below in issue report
                # investigate f-string formatting here
                msg = (
                    'Warning: Metadata file keyword TYPE provided as '
                    '{self.metadata_types[0]}'
                )
                self.store_validation_issue('warn', 'format', msg)
        else:
            msg = 'Error: Metadata file TYPE row malformed, missing TYPE'
            self.store_validation_issue('error', 'format', msg)
        return valid

    def validate_type_annotations(self):
        """Check metadata second row contains only 'group' or 'numeric'.
        :return: boolean   True if all type annotations are valid, otherwise False
        """
        valid = False
        invalid_types = []
        # skipping the TYPE keyword, iterate through the types
        # collecting invalid type annotations in list annots
        for t in set(list(self.file.columns)[1:]):
            if t[1] not in self.annotation_type:
                # if the value is a blank space, store a higher visibility
                # string for error reporting
                if not t:
                    invalid_types.append('<empty value>')
                else:
                    invalid_types.append(t)
        if invalid_types:
            msg = 'Error: TYPE declarations should be group or numeric'
            self.store_validation_issue('error', 'format', msg, invalid_types)
        else:
            valid = True
        return valid

    def validate_against_header_count(self):
        """Metadata header and type counts should match.
        :return: boolean   True if header and type counts match, otherwise False
        """
        valid = False
        len_headers = len(self.file.columns.labels[0])
        len_annot_type = len(self.file.columns.labels[1])
        if not len_headers == len_annot_type:
            msg = (
                f'Error: {len_annot_type} TYPE declarations '
                f'for {len_headers} column headers'
            )
            self.store_validation_issue('error', 'format', msg)
        else:
            valid = True
        return valid

    def validate_format(self):
        """Check all metadata file format criteria for file validity
        """
        return (
            self.validate_header_keyword()
            or self.validate_type_keyword()
            or self.validate_type_annotations()
            or self.validate_unique_header()
            or self.validate_against_header_count()
        )
