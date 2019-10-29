"""Module for text, CSV, and TSV file types used in SCP ingest

DESCRIPTION
Module provides extract capabilities for text, CSV, and TSV file types

"""
import csv
import mimetypes
import os
import re
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
from dataclasses import dataclass
from mypy_extensions import TypedDict

import pandas as pd
from google.cloud import storage


@dataclass
class DataArray(TypedDict):
    MAX_ENTREIS = 100_000

    def __init__(
        self,
        name: str,
        cluster_name: str,
        array_type: str,
        values: List,
        linear_data_type: str,
        linear_data_id: str,
        study_id: str,
        study_file_id: str,
        array_index: int = 0,
        subsample_threshold: int = None,
        subsample_annotation: str = None,
    ):
        self.name = name
        self.cluster_name = cluster_name
        self.array_type = array_type
        self.array_index = array_index
        self.values = values
        self.subsample_threshold = subsample_threshold
        self.subsample_annotation = subsample_annotation
        self.linear_data_type = linear_data_type
        self.linear_data_id = linear_data_id
        self.study_id = study_id
        self.study_file_id = study_file_id

        # TODO: Add logic for when len(self.values) > self.MAX_ENTREIS
        # get_dataArray(self):
        #     if len(self.values) > self.MAX_ENTREIS:


class IngestFiles:
    def __init__(self, file_path, allowed_file_types, *, open_as=None):
        self.file_path = file_path
        # File is remote (in GCS bucket) when running via PAPI,
        # and typically local when developing
        self.is_remote_file = file_path[:5] == "gs://"

        self.verify_file_exists(file_path)

        self.allowed_file_types = allowed_file_types
        self.file_type, self.file, self.file_handle = self.open_file(file_path, open_as)
        # Keeps tracks of lines parsed
        self.amount_of_lines = 0

    def download_from_bucket(self, file_path):
        """Downloads file from Google Cloud Storage bucket"""
        blob = self.bucket.blob(self.source)
        destination = "/tmp/" + self.source.replace("/", "%2f")
        blob.download_to_filename(destination)
        print(f"{file_path} downloaded to {destination}.")
        return destination

    def set_gcs_attrs(self, file_path):
        """Sets instance attributes related to Google Cloud Storage"""
        self.storage_client = storage.Client()
        path_segments = file_path[5:].split("/")
        self.bucket_name = path_segments[0]
        self.bucket = self.storage_client.get_bucket(self.bucket_name)
        self.source = "/".join(path_segments[1:])

    def verify_file_exists(self, file_path):
        """Determines if file can be found, throws error if not"""
        if self.is_remote_file:
            # File is in GCS bucket
            self.set_gcs_attrs(file_path)
            source_blob = storage.Blob(bucket=self.bucket, name=self.source)
            if not source_blob.exists(self.storage_client):
                raise OSError(f'Remote file "{file_path}" not found')
        else:
            # File is local
            if not os.path.exists(file_path):
                raise OSError(f'File "{file_path}" not found')

    def resolve_path(self, file_path):
        """Localizes object if given a GS URL, returns open Python file object

        Args:
            file_path: Path to a local file, or a Google Cloud Storage URL

        Returns:
            Open file object
        """
        if self.is_remote_file:
            file_path = self.download_from_bucket(file_path)

        # Remove BOM with encoding ='utf - 8 - sig'
        return open(file_path, encoding="utf-8-sig")

    def reset_file(self, start_point, open_as=None):
        """Restart file reader at point that's equal to start_point.
        Method is used in cases where a file may need to be read multiple times"""

        self.file_type, self.file, self.file_handle = self.open_file(
            self.file_path, start_point=start_point, open_as=open_as
        )

    def open_file(self, file_path, open_as=None, start_point: int = 0):
        """ Opens txt, csv, or tsv formatted files"""
        open_file = self.resolve_path(file_path)
        if start_point != 0:
            for i in range(start_point):
                open_file.readline()
        file_connections = {
            "text/csv": self.open_csv(open_file),
            "text/plain": open_file,
            "application/json": open_file,
            "text/tab-separated-values": self.open_tsv(open_file),
            "dataframe": self.open_pandas,
        }
        # Check file type
        file_type = self.get_file_type(file_path)[0]
        # See if file type is allowed
        if file_type in self.allowed_file_types:
            # Return file object and type
            if open_as is None:
                return file_type, file_connections.get(file_type), open_file
            else:
                return (
                    file_type,
                    file_connections.get("dataframe")(open_file, file_path),
                    open_file,
                )
        else:
            raise ValueError(
                f"Unsupported file format. Allowed file types are: {' '.join(self.allowed_file_type)}"
            )

    # Inherited function
    def extract(self):
        """ Calls extract function for txt, csv, or tsv formatted files to
            retrieve all contents from file.
        """

        file_type_extract_fns = {
            "text/csv": self.extract_csv_or_tsv,
            "text/plain": self.extract_txt,
            "text/tab-separated-values": self.extract_csv_or_tsv,
        }
        return file_type_extract_fns.get(self.file_type)()

    def split_line(self, line):
        """Splits lines on file format-appropriate delimiters"""

        return re.findall(r"[^,\t]+", line)

    def get_file_type(self, file_path):
        """Returns file type"""
        return mimetypes.guess_type(file_path)

    def open_pandas(self, opened_file, file_path):
        """Opens file as a dataframe """
        opened_file.readline()
        meta_data = opened_file.readline()
        if meta_data.find("\t") != -1:
            return pd.read_csv(
                file_path, sep="\t", header=[0, 1], quoting=csv.QUOTE_NONE
            )
        elif meta_data.find(",") != -1:
            return pd.read_csv(
                file_path, sep=",", header=[0, 1], quoting=csv.QUOTE_NONE
            )
        else:
            raise ValueError("File must be tab or comma delimited")

    def merge_df(self, first_df, second_df):
        """ Does an inner join on a dataframe """
        self.file = pd.merge(second_df, first_df, on=[("NAME", "TYPE")])

    def open_csv(self, opened_file_object):
        """Opens csv file"""
        csv.register_dialect(
            "csvDialect", delimiter=",", quoting=csv.QUOTE_ALL, skipinitialspace=True
        )
        return csv.reader(opened_file_object, dialect="csvDialect")

    def open_tsv(self, opened_file_object):
        """Opens tsv file"""
        csv.register_dialect(
            "tsvDialect", delimiter="\t", quoting=csv.QUOTE_ALL, skipinitialspace=True
        )
        return csv.reader(opened_file_object, dialect="tsvDialect")

    def extract_csv_or_tsv(self):
        """Extracts all rows from a csv or tsv file"""
        while True:
            try:
                row = next(self.file)
                self.amount_of_lines += 1
                return row
            except StopIteration:
                break

    def extract_txt(self):
        """Extracts all lines from txt files

        Returns:
                next_row_revised : List[str]
                    A single row from a txt file.
        """
        while True:
            next_row = self.file.readline()
            if not next_row:
                break
            self.amount_of_lines += 1
            # Create array with no new line, commas, or tab characters
            next_row_revised = self.split_line(
                next_row.replace("\n", "").replace('"', "")
            )
            return next_row_revised

    def get_next_line(self, *, increase_line_count=True, split_line=True):
        """Returns a single line of txt, csv or tsv files"""

        next_row = next(self.file)
        # Increase counter for line extracted
        if increase_line_count:
            self.amount_of_lines += 1
        elif self.file_type in ("text/csv", "text/tab-separated-values"):
            # CSV and TSV files returns next lines as array split on tabs/commas
            return next_row
        elif self.file_type == "text/plain":
            # Create array with no new line, commas, or tab characters
            next_row_revised = next_row.replace("\n", "")
            if split_line:
                return self.split_line(next_row_revised)
            else:
                return next_row_revised

    # This function will be deleted and in replaced in annotations.py
    def determine_coordinates_and_cell_names(self):
        """Finds column names for coordinates, annotations, and cell names"""
        self.coordinates_and_cell_names = [
            annot[0]
            for annot in self.file.columns
            if annot[0].lower() in ('z', 'y', 'x', 'name')
        ]
        # annotation column names
        self.columns = [
            annot
            for annot in self.file.columns
            if annot[0].lower() not in ('z', 'y', 'x', 'name')
        ]

    # This function will be deleted and in replaced in annotations.py
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
