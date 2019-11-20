"""Base class for annotation files (such as cluster and metadata files)

DESCRIPTION
Class defines common functions and validation methods for annotation files

PREREQUISITES
Must have python 3.6 or higher.
"""

import abc
import pandas as pd  # NOqa: F821
from ingest_files import IngestFiles


class Annotations(IngestFiles):
    __metaclass__ = abc.ABCMeta

    def __init__(
        self, file_path, allowed_file_types, study_id=None, study_file_id=None
    ):
        IngestFiles.__init__(
            self, file_path, allowed_file_types, open_as='dataframe', header=[0, 1]
        )
        self.headers = self.file.columns.get_level_values(0)
        self.annot_types = self.file.columns.get_level_values(1)
        self.study_id = study_id
        self.study_file_id = study_file_id

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

    def preproccess(self):
        """Ensures that:
            - Numeric columns are rounded to 3 decimals points
            - Group annotations are strings
            - 'NAME' in first header row is capitalized
            - 'TYPE' in second header row is capitalized
        """
        headers = self.file.columns.get_level_values(0)
        annot_types = self.file.columns.get_level_values(1)
        # Lowercase second level. Example: NUMeric -> numeric
        self.file.rename(
            columns=lambda col_name: col_name.lower(), level=1, inplace=True
        )
        name = list(headers)[0]
        type = list(annot_types)[0].lower()
        # Uppercase NAME and TYPE
        self.file.rename(columns={name: name.upper(), type: type.upper()}, inplace=True)
        # Make sure group annotations are treated as strings
        group_columns = self.file.xs(
            "group", axis=1, level=1, drop_level=False
        ).columns.tolist()
        self.file[group_columns] = self.file[group_columns].astype(str)
        # Find numeric columns and round to 3 decimals places and are floats
        numeric_columns = self.file.xs(
            "numeric", axis=1, level=1, drop_level=False
        ).columns.tolist()
        # TODO perform replace
        self.file[numeric_columns] = self.file[numeric_columns].round(3).astype(int)

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
        """Check all metadata file format criteria for file validity
        """
        self.is_valid_file = all(
            [
                self.validate_header_keyword(),
                self.validate_type_keyword(),
                self.validate_type_annotations(),
                self.validate_unique_header(),
                self.validate_against_header_count(),
            ]
        )
        return self.is_valid_file
