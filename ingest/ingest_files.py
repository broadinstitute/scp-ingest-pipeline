"""Module for text, CSV, and TSV file types used in SCP ingest

DESCRIPTION
Module provides extract capabilities for text, CSV, and TSV file types

"""
import csv
import mimetypes
import os
import re
from itertools import islice


class IngestFiles:
    def __init__(self, file_path, allowed_file_types, is_MTX=False):
        if not os.path.exists(file_path):
            raise IOError(f"File '{file_path}' not found")
        self.allowed_file_types = allowed_file_types
        self.file_type, self.file = self.open_file(file_path)
        # Keeps tracks of lines parsed
        self.amount_of_lines = 0
        self.is_MTX = is_MTX

    def open_file(self, file_path):
        """ Opens txt, csv, or tsv formatted files"""
        open_file = self.resolve_path(file_path)
        file_connections = {
            # Remove BOM with encoding='utf-8-sig'
            'text/csv': self.open_csv(open_file),
            'text/plain': open_file,
            'text/tab-separated-values': self.open_tsv(open_file),
        }
        # Check file type
        file_type = self.get_file_type(file_path)[0]
        # See if file type is allowed
        if file_type in self.allowed_file_types:
            # Return file object and type
            return file_type, file_connections.get(file_type)
        else:
            raise ValueError(f"Unsupported file format. Allowed file types are: {' '.join(self.allowed_file_type)}")

    def resolve_path(self, file_path):
        """Localizes object if given a GS URL, returns open Python file object

        Args:
            file_path: Path to a local file, or a Google Cloud Storage URL

        Returns:
            Open file object
        """
        if file_path[:5] == 'gs://':
            split_path = file_path[5:].split('/')
            bucket_name = split_path[0]
            source = '/'.join(split_path[1:])
            destination = '/tmp/' + source.replace('/', '%2f')
            download_blob(bucket_name, source, destination)
            file_path = destination

        return open(file_path, encoding='utf-8-sig')

    # Inherited function
    def extract(self):
        """ Calls extract function for txt, csv, or tsv formatted files to
            retrieve all contents from file.
        """

        file_type_extract_fns = {
            'text/csv': self.extract_csv_or_tsv,
            'text/plain': self.extract_txt,
            'text/tab-separated-values': self.extract_csv_or_tsv,
        }
        return file_type_extract_fns.get(self.file_type)()

    def split_line(self, line):
        """Splits lines on file format-appropriate delimiters"""

        if self.is_MTX:
            return re.findall(r'[^,\t]+', line)
        else:
            return re.findall(r'[^,\s\t]+', line)

    def get_file_type(self, file_path):
        """Returns file type"""
        return mimetypes.guess_type(file_path)

    def open_csv(self, opened_file_object):
        """Opens csv file"""
        csv.register_dialect('csvDialect',
                             delimiter=',',
                             quoting=csv.QUOTE_ALL,
                             skipinitialspace=True)
        return csv.reader(opened_file_object, dialect='csvDialect')

    def open_tsv(self, opened_file_object):
        """Opens tsv file"""
        csv.register_dialect('tsvDialect',
                             delimiter='\t',
                             quoting=csv.QUOTE_ALL,
                             skipinitialspace=True)
        return csv.reader(opened_file_object, dialect='tsvDialect')

    def extract_csv_or_tsv(self):
        """Extracts all rows from a csv or tsv file"""
        while(True):
            try:
                row = next(self.file)
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
            next_row_revised = self.split_line(next_row.replace('\n', ''))
            return next_row_revised

    def get_next_line(self, *, increase_line_count=True, split_line=True):
        """Returns a single line of txt, csv or tsv files"""

        next_row = next(self.file)
        # Increase counter for line extracted
        if increase_line_count:
            self.amount_of_lines += 1
        elif self.file_type in ('text/csv', 'text/tab-separated-values'):
            # CSV and TSV files returns next lines as array split on tabs/commas
            return next_row
        elif self.file_type == 'text/plain':
            # Create array with no new line, commas, or tab characters
            next_row_revised = next_row.replace('\n', '')
            if split_line:
                return self.split_line(next_row_revised)
            else:
                return next_row_revised
