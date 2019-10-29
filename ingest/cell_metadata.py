"""Module for ingesting cell metadata files
DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.
PREREQUISITES
Must have python 3.6 or higher.
"""
from collections import defaultdict
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
from dataclasses import dataclass
from mypy_extensions import TypedDict
import ntpath

from ingest_files import DataArray
from annotations import Annotations

try:
    # Used when importing internally and in tests
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles

DOCUMENT_LIMIT_BYTES = 1_048_576


class CellMetadata(Annotations):
    ALLOWED_FILE_TYPES = ['text/csv', 'text/plain', 'text/tab-separated-values']

    def __init__(
        self, file_path: str, study_id: str, study_file_id: str, *args, **kwargs
    ):

        IngestFiles.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, open_as='dataframe'
        )
        self.preproccess()
        self.file_path = file_path
        self.headers = self.file.columns.get_level_values(0)
        self.annot_types = self.file.columns.get_level_values(1)
        self.cell_names = []
        self.study_id = study_id
        self.study_file_id = study_file_id
        # lambda below initializes new key with nested dictionary as value and avoids KeyError
        self.issues = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.ontology = defaultdict(lambda: defaultdict(list))
        self.cells = []
        self.is_valid_file = self.validate_format()

    # This model pertains to columns from cell metadata files
    @dataclass
    class Model(TypedDict):
        # value from column header
        name: str
        annotation_type: str
        # unique values from "group" type annotations
        values: List
        study_file_id: str
        study_id: str

    def transform(self):
        """ Transform data from cell metadata files into data model"""
        # first column is cell names, therefore skip
        for annot_header in self.file.columns[1:]:
            annot_name = annot_header[0]
            annot_type = annot_header[1]
            yield self.Model(
                {
                    'name': annot_name,
                    'study_id': self.study_id,
                    # unique values from "group" type annotations else []
                    'values': list(self.file[annot_header].unique())
                    if annot_type == 'group'
                    else [],
                    'annotation_type': annot_type,
                    'study_file_id': self.study_file_id,
                }
            )

    def set_data_array(self, annot_header: str, linear_data_id: str):
        data_array_attrs = locals()
        annot_name = annot_header[0]
        head, tail = ntpath.split(self.file_path)
        base_data_array_model = {
            'cluster_name': tail or ntpath.basename(head),
            'value': list(self.file[annot_name]),
            'study_file_id': self.study_file_id,
            'study_id': self.study_id,
        }
        # This is an array (or group of arrays) of every cell
        if annot_name.lower() == 'name':
            base_data_array_model.update(
                {
                    'name': 'All Cells',
                    'array_type': 'cells',
                    'linear_data_type': 'Study',
                }
            )
        # data from cell metadata file that correspond to a column of data
        else:
            base_data_array_model.update(
                {
                    'name': 'annot_name',
                    'array_type': 'annotations',
                    'linear_data_type': 'CellMetadatum',
                }
            )
        return DataArray({**locals(), **base_data_array_model})

    def yield_by_row(self) -> None:
        """ Yield row from cell metadata file"""
        for row in self.file.itertuples(index=False):
            dict_row = row._asdict()
            yield dict_row.values()

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
        if self.headers[0].upper() == 'NAME':
            valid = True
            if self.headers[0] != 'NAME':
                msg = f'Metadata file keyword "NAME" provided as ' f"{self.headers[0]}"
                self.store_validation_issue('warn', 'format', msg)
        else:
            msg = 'Malformed metadata file header row, missing NAME. (Case Sensitive)'
            self.store_validation_issue('error', 'format', msg)
        return valid

    def validate_unique_header(self):
        """Check all metadata header names are unique and not empty.
        :return: boolean   True if valid, False otherwise
        """
        valid = False
        unique_headers = set(self.headers)
        if len(unique_headers) == len(self.headers):
            valid = True
        else:
            seen_headers = set()
            duplicate_headers = set()
            for x in self.headers:
                if x in seen_headers or seen_headers.add(x):
                    duplicate_headers.add(x)
            msg = (
                f'Duplicated metadata header names are not allowed: {duplicate_headers}'
            )
            self.store_validation_issue('error', 'format', msg)
            valid = False
        if any('Unnamed' in s for s in list(unique_headers)):
            msg = 'Headers cannot contain empty values'
            self.store_validation_issue('error', 'format', msg)
            valid = False
        return valid

    def validate_type_keyword(self):
        """Check metadata second row starts with TYPE (case-insensitive).
        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.annot_types[0].upper() == 'TYPE':
            valid = True
            if self.annot_types[0] != 'TYPE':
                msg = f'Metadata file keyword "TYPE" provided as {self.annot_types[0]}'
                self.store_validation_issue('warn', 'format', msg)
        else:
            msg = 'Malformed metadata TYPE row, missing TYPE. (Case Sensitive)'
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
        for t in self.annot_types[1:]:
            if t.lower() not in ('group', 'numeric'):
                # if the value is a blank space, store a higher visibility
                # string for error reporting
                if 'Unnamed' in t:
                    invalid_types.append('<empty value>')
                # Duplicated metadata header name causes type annotation issue.
                # Side effect of Pandas adding a suffix to uniquefy the header.
                # These invalid annotations should not be included in invalid
                # type annotation count. This exception may cause miscount of
                # type annot errors if user-supplied annotation has period.
                elif '.' in t:
                    pass
                else:
                    invalid_types.append(t)
        if invalid_types:
            msg = 'TYPE row annotations should be "group" or "numeric"'
            self.store_validation_issue('error', 'format', msg, invalid_types)
        else:
            valid = True
        return valid

    def validate_against_header_count(self):
        """Metadata header and type counts should match.
        :return: boolean   True if header and type counts match, otherwise False
        """
        valid = False
        len_headers = len(
            [header for header in self.headers if 'Unnamed' not in header]
        )
        len_annot_type = len(
            [
                annot_type
                for annot_type in self.annot_types
                if 'Unnamed' not in annot_type
            ]
        )
        if not len_headers == len_annot_type:
            msg = (
                f'Header mismatch: {len_annot_type} TYPE declarations '
                f'for {len_headers} column headers'
            )
            self.store_validation_issue('error', 'format', msg)
        else:
            valid = True
        return valid

    def validate_format(self):
        """Check all metadata file format criteria for file validity
        """
        return all(
            [
                self.validate_header_keyword(),
                self.validate_type_keyword(),
                self.validate_type_annotations(),
                self.validate_unique_header(),
                self.validate_against_header_count(),
            ]
        )
