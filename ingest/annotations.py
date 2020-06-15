"""Base class for annotation files (such as cluster and metadata files)

DESCRIPTION
Class defines common functions and validation methods for annotation files

PREREQUISITES
Must have python 3.6 or higher.
"""

import abc
from collections import defaultdict

import pandas as pd
from bson.objectid import ObjectId

try:
    # Used when importing internally and in tests
    from ingest_files import IngestFiles
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from .ingest_files import IngestFiles


class Annotations(IngestFiles):
    __metaclass__ = abc.ABCMeta

    def __init__(
        self, file_path, allowed_file_types, study_id=None, study_file_id=None
    ):
        IngestFiles.__init__(self, file_path, allowed_file_types)
        if study_id is not None:
            self.study_id = ObjectId(study_id)
        if study_file_id is not None:
            self.study_file_id = ObjectId(study_file_id)
        # lambda below initializes new key with nested dictionary as value and avoids KeyError
        self.issues = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        csv_file, self.file_handle = self.open_file(self.file_path)
        # Remove white spaces, quotes (only lowercase annot_types)
        self.headers = [header.strip().strip('\"') for header in next(csv_file)]
        self.annot_types = [type.strip().strip('\"').lower() for type in next(csv_file)]

    def reset_file(self):
        self.preprocess()

    @abc.abstractmethod
    def transform(self):
        """Returns data model"""

    @abc.abstractmethod
    def set_data_array(self):
        """Sets DataArray"""

    def determine_coordinates_and_cell_names(self):
        """Finds column names for coordinates, annotations, and cell names"""
        self.coordinates_and_cell_headers = [
            annot[0]
            for annot in self.file.columns
            if annot[0].lower() in ('z', 'y', 'x', 'name')
        ]
        # annotation column names
        self.annot_column_headers = [
            annot
            for annot in self.file.columns
            if annot[0].lower() not in ('z', 'y', 'x', 'name')
        ]

    def merge_df(self, first_df, second_df):
        """ Does an inner join on a dataframe """
        self.file = pd.merge(second_df, first_df, on=[("NAME", "TYPE")])

    def preprocess(self):
        """Ensures that:
            - 'NAME' in first header row is capitalized
            - 'TYPE' in second header row is capitalized
        """

        # Uppercase NAME and TYPE
        self.headers[0] = self.headers[0].upper()
        self.annot_types[0] = self.annot_types[0].upper()
        self.create_data_frame()

    def create_data_frame(self):
        """
        - Create dataframe with proper dtypes to ensure:
            - Labels are treated as strings (objects)
            - Numeric annotations are treated as float32
            - Numeric columns are rounded to 3 decimals points
        """
        columns = []
        dtypes = {}
        for annotation, annot_type in zip(self.headers, self.annot_types):
            dtypes[annotation] = 'object' if annot_type != 'numeric' else 'float32'
            columns.append((annotation, annot_type))
            # multiIndex = pd.MultiIndex.from_tuples(index)
        index = pd.MultiIndex.from_tuples(columns)
        self.file = self.open_file(
            self.file_path, open_as='dataframe', dtype=dtypes, names=index, skiprows=2
        )[0]
        # dtype of object allows mixed dtypes in columns, including numeric dtypes
        # coerce group annotations that pandas detects as non-object types to type string
        for annotation, annot_type in columns:
            if (
                annot_type == 'group'
                and self.file[annotation].dtypes[annot_type] != 'O'
            ):
                self.file[annotation] = self.file[annotation].astype(str)
        if 'numeric' in self.annot_types:
            numeric_columns = self.file.xs(
                "numeric", axis=1, level=1, drop_level=False
            ).columns.tolist()
            try:
                # Round numeric columns to 3 decimal places
                self.file[numeric_columns] = (
                    self.file[numeric_columns].round(3).astype(float)
                )
            except Exception as e:
                self.error_logger.error(
                    "There are non-numeric values in numeric columns",
                    extra=self.extra_log_params,
                )
                self.error_logger.error(e, extra=self.extra_log_params)

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
        """Check header row starts with NAME (case-insensitive).

            return: boolean   True if valid, False otherwise
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
        """Check all header names are unique and not empty.
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
        """Check second row starts with TYPE (case-insensitive).
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
        """Ensure second row contains only 'group' or 'numeric'.
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
        """Header and type counts should match.
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
        """Check common format criteria for annotation files
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
