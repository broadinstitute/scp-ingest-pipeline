"""Module for text, CSV, and TSV file types used in SCP ingest

DESCRIPTION
Module provides extract capabilities for text, CSV, and TSV file types
"""
import copy
import csv
import gzip
import logging
import mimetypes
import os
import re
from dataclasses import dataclass
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401

import pandas as pd  # NOqa: F821
from google.cloud import storage

# from google.cloud.logging.resource import Resource
# import google.cloud.logging

try:
    from monitor import setup_logger
except ImportError:
    from .monitor import setup_logger


@dataclass
class DataArray:
    MAX_ENTRIES = 100_000
    COLLECTION_NAME = "data_arrays"
    errors_logger = setup_logger(
        __name__ + "_errors", "errors.txt", level=logging.ERROR
    )

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
        # special case to override linear_data_id in case of 'Study' linear_data_type
        if self.linear_data_type == "Study":
            self.linear_data_id = self.study_id

    def get_data_array(self):
        if len(self.values) > self.MAX_ENTRIES:
            values = self.values
            for idx, i in enumerate(range(0, len(self.values), self.MAX_ENTRIES)):
                self.values = values[i : i + self.MAX_ENTRIES]
                self.array_index = idx
                yield copy.copy(self.__dict__)
        else:
            yield self.__dict__


class IngestFiles:
    # General logger for class
    info_logger = setup_logger(__name__, "info.txt")
    error_logger = setup_logger(__name__ + "_errors", "errors.txt", level=logging.ERROR)

    def __init__(self, file_path, allowed_file_types):
        self.file_path = file_path
        # File is remote (in GCS bucket) when running via PAPI,
        # and typically local when developing
        self.is_remote_file = IngestFiles.is_remote_file(file_path)
        self.is_gzip_file = self.get_file_type(file_path)[1] == 'gzip'

        self.verify_file_exists(file_path)
        # Allowed files for a given file type (expression file, cluster files, etc.)
        self.allowed_file_types = allowed_file_types
        # Keeps tracks of lines parsed
        self.amount_of_lines = 0

    @staticmethod
    def is_remote_file(file_path):
        return file_path[:5] == "gs://"

    def download_from_bucket(self, file_path):
        """Downloads file from Google Cloud Storage bucket"""
        blob = self.bucket.blob(self.source)
        destination = "/tmp/" + self.source.replace("/", "%2f")
        blob.download_to_filename(destination)
        self.info_logger.info(
            f"{file_path} downloaded to {destination}.",
            extra={"study_id": None, "duration": None},
        )
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
                self.error_logger.error(
                    f'Remote file "{file_path}" not found',
                    extra={"study_id": None, "duration": None},
                )
                raise OSError(f'Remote file "{file_path}" not found')
        else:
            # File is local
            if not os.path.exists(file_path):
                self.error_logger.error(
                    f'File "{file_path}" not found',
                    extra={"study_id": None, "duration": None},
                )
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
            self.local_file_path = file_path
        # Remove BOM with encoding ='utf - 8 - sig'
        if self.is_gzip_file:
            open_file = gzip.open(file_path, "rt", encoding="utf-8-sig")
        else:
            open_file = open(file_path, 'rt', encoding="utf-8-sig")
        return open_file, file_path

    def reset_file(self, file_path, start_point, open_as=None):
        """Restart file reader at point that's equal to start_point.
        Method is used in cases where a file may need to be read multiple times"""

        return self.open_file(file_path, start_point=start_point, open_as=open_as)

    @staticmethod
    def delocalize_file(
        study_file_id, study_id, file_path, file_to_delocalize, bucket_destination
    ):
        """Writes local file to Google bucket
        Args:
            file_path: path of an ingest file (MUST BE  GS url)
            file_to_delocalize: name of local file to delocalize (ie. errors.txt)
            bucket_destination: path to google bucket (ie. parse_logs/{study_file_id}/errors.txt)

        """
        info_logger = setup_logger(__name__, "info.txt")
        error_logger = setup_logger(
            __name__ + "_errors", "errors.txt", level=logging.ERROR
        )
        extra_log_params = {"study_id": study_id, "duration": None}
        if IngestFiles.is_remote_file(file_path):
            try:
                path_segments = file_path[5:].split("/")
                bucket_name = path_segments[0]
                storage_client = storage.Client()
                bucket = storage_client.get_bucket(bucket_name)
                blob = bucket.blob(bucket_destination)
                blob.upload_from_filename(file_to_delocalize)
                info_logger.info(
                    f"File {file_to_delocalize} uploaded to {bucket_destination}.",
                    extra=extra_log_params,
                )
            except Exception as e:
                error_logger.error(
                    f"File {file_to_delocalize} not uploaded to {bucket_destination}.",
                    extra=extra_log_params,
                )
                error_logger.error(e, extra=extra_log_params)
        else:
            error_logger.error(
                "Cannot push to bucket. File is not remote", extra=extra_log_params
            )

    def open_file(self, file_path, open_as=None, start_point: int = 0, **kwargs):
        """ A wrapper function for opening txt (txt is expected to be tsv or csv), csv, or tsv formatted files"""
        open_file, file_path = self.resolve_path(file_path)
        file_connections = {
            "text/csv": self.open_csv,
            "text/plain": self.open_txt,
            "text/tab-separated-values": self.open_tsv,
            "dataframe": self.open_pandas,
        }

        if start_point != 0:
            self.amount_of_lines = 0
            for i in range(start_point):
                open_file.readline()
        # See if file type is allowed
        file_type = self.get_file_type(file_path)[0]
        self.info_logger.info(
            f"opening {file_path} as: {file_type}",
            extra={"study_id": None, "duration": None},
        )
        if file_type in self.allowed_file_types:
            # Return file object and type
            if file_type == "application/json":
                return open_file
            elif open_as is None:
                if file_type == "text/plain":
                    return (
                        file_connections.get(file_type)(open_file, file_type, **kwargs),
                        open_file,
                    )
                else:
                    return (
                        file_connections.get(file_type)(open_file, **kwargs),
                        open_file,
                    )
            else:
                return (
                    file_connections.get(open_as)(
                        file_path, file_type, open_file_object=open_file, **kwargs
                    ),
                    open_file,
                )
        else:
            raise ValueError(
                f"Unsupported file format. Allowed file types are: {' '.join(self.allowed_file_types)}"
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

    def open_txt(self, open_file_object, file_type, **kwargs):
        """Method for opening txt files that are expected be tab
        or comma delimited"""
        if file_type == "text/tab-separated-values":
            delimiter = "\t"
        elif file_type == "text/csv":
            delimiter = ","
        else:
            delimiter = None
        # Determine if file is tsv or csv
        # reading single line of file instead of first 1024 bytes to avoid issues with delimiter detection
        # reference: https://stackoverflow.com/questions/35756682/getting-csv-sniffer-to-work-with-quoted-values
        try:
            csv_dialect = csv.Sniffer().sniff(
                open_file_object.readline(), delimiters=delimiter
            )
            csv_dialect.skipinitialspace = True
            open_file_object.seek(0)
            return csv.reader(open_file_object, csv_dialect)
        except Exception:
            raise ValueError(
                f'Could not determine delimiter. Please save file with appropriate suffix (.tsv or .csv) and try again.'
            )

    def open_pandas(self, file_path, file_type, **kwargs):
        """Opens file as a dataframe """
        open_file_object = kwargs.pop('open_file_object')
        if file_type in self.allowed_file_types:
            # Determine delimiter based on file type
            if file_type == "text/tab-separated-values":
                delimiter = "\t"
            elif file_type == "text/csv":
                delimiter = ","
            else:
                delimiter = None
            dialect = csv.Sniffer().sniff(
                open_file_object.readline(), delimiters=delimiter
            )
            dialect.skipinitialspace = True
            open_file_object.seek(0)
            return pd.read_csv(file_path, dialect=dialect, low_memory=False, **kwargs)
        else:
            raise ValueError("File must be tab or comma delimited")

    def open_csv(self, opened_file_object, **kwargs):
        """Opens csv file"""
        csv.register_dialect(
            "csvDialect",
            delimiter=",",
            quotechar='"',
            skipinitialspace=True,
            escapechar='\\',
        )
        return csv.reader(opened_file_object, dialect="csvDialect")

    def open_tsv(self, opened_file_object, **kwargs):
        """Opens tsv file"""
        csv.register_dialect("tsvDialect", delimiter="\t", skipinitialspace=True)
        return csv.reader(opened_file_object, dialect="tsvDialect")

    def extract_csv_or_tsv(self, file):
        """Extracts all rows from a csv or tsv file"""
        while True:
            try:
                row = next(file)
                self.amount_of_lines += 1
                return row
            except StopIteration:
                break
