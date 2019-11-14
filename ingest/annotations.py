"""Class for annotation files (such as cluster and metadata files)

DESCRIPTION
Class defines common functions needed for annotation type files

PREREQUISITES
Must have python 3.6 or higher.

"""

import abc

from ingest_files import IngestFiles


class Annotations(IngestFiles):
    __metaclass__ = abc.ABCMeta

    def __init__(self, file_path, allowed_file_types):
        IngestFiles.__init__(
            self, file_path, allowed_file_types, open_as='dataframe', header=[0, 1]
        )

    @abc.abstractmethod
    def transform(self):
        """Returns data model"""

    # This will end up being a class method
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
        self.file[numeric_columns] = self.file[numeric_columns].round(3).astype(float)
