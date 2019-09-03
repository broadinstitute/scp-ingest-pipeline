"""Module for ingesting cell metadata files
DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.
PREREQUISITES
Must have python 3.6 or higher.
"""

from collections import defaultdict
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
import warnings

from ingest_files import IngestFiles


class CellMetadata(IngestFiles):
    ALLOWED_FILE_TYPES = ["text/csv", "text/plain", "text/tab-separated-values"]
    COLLECTION_NAME = "cell_metadata"
    SUBCOLLECTION_NAME = "data"

    def __init__(self, file_path, file_id: str, study_accession: str, *args, **kwargs):

        IngestFiles.__init__(
            self, file_path, self.ALLOWED_FILE_TYPES, open_as="dataframe"
        )
        # unique values for group-based annotations
        self.unique_values = []
        self.cell_names = []
        self.annotation_type = ["group", "numeric"]
        self.study_accession = study_accession
        self.file_id = file_id
        self.errors = defaultdict(list)
        self.validate_format()
        self.ontology = defaultdict(lambda: defaultdict(set))
        self.type = defaultdict(list)
        self.cells = []
        self.preproccess()

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
                {
                    "name": col_name.lower()
                    if col_name.lower() in ("x", "y", "z")
                    else col_name,
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

    def validate_header_keyword(self):
        """Check metadata header row starts with NAME (case-insensitive).

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.file.columns[0][0].upper() == "NAME":
            valid = True
            if self.file.columns[0] != "NAME":
                # ToDO - capture warning below in error report
                warnings.warn(
                    f'Warning: metadata file keyword "NAME" provided as '
                    f"{self.file.columns[0][0]}",
                    UserWarning,
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
        if len(set(self.file.columns.labels[0])) == len(self.file.columns.labels[0]):
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
        if self.file.columns[0][1].upper() == "TYPE":
            valid = True
            if self.file.columns[0][1] != "TYPE":
                # ToDO - capture warning below in error report
                # investigate f-string formatting here
                warnings.warn(
                    f'Warning: metadata file keyword "TYPE" provided as {self.file.columns[0][1]}',
                    UserWarning,
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
        for t in set(list(self.file.columns)[1:]):
            if t[1] not in self.annotation_type:
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
        if not len(self.file.columns.levels[0]) == len(self.file.columns.levels[1]):
            self.errors["format"].append(
                f"Error: {len(self.file.columns.levels[1])} TYPE declarations "
                f"for {len(self.file.columns.levels[0])} column headers"
            )
        else:
            valid = True
        return valid

    def validate_empty_header(self):
        """Check for empty values in header"""
        valid = False
        if len(self.file.filter(like="Unnamed:").columns.values) > 0:
            self.errors["format"].append("Error: Empty value in header.")
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
        self.validate_empty_header()
        if self.errors["format"]:
            valid = False
        else:
            valid = True
        return valid
