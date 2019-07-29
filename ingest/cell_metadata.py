"""Module for ingesting cell metadata files

DESCRIPTION
Module provides extract and transform functions for cell metadata files.
Text, CSV, and TSV files are supported.

PREREQUISITES
Must have python 3.6 or higher.
"""
import copy
import json
import os
from typing import Dict, Generator, List, Tuple, Union

from ingest_files import IngestFiles


class CellMetadata(IngestFiles):
    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']

    def __init__(self, file_path, file_id: str = None, study_accession: str = None):

        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.header = self.get_next_line(
            increase_line_count=False)
        self.metadata_types = self.get_next_line(
            increase_line_count=False)
        # unique values for group-based annotations
        self.uniqueValues = []
        self.cell_names = []
        self.top_level_doc = self.create_documents(
            file_path, file_id, study_accession)
        self.data_subcollection = self.create_subdocuments()

    def transform(self, row: List[str]) -> None:
        """ Add data from cell metadata files into data model"""
        for idx, column in enumerate(row):
            if idx != 0:
                # if annotation is numeric convert from string to float
                if self.metadata_types[idx].lower() == 'numeric':
                    column = round(float(column), 3)
                elif self.metadata_types[idx].lower() == 'group':
                    # Check for unique values
                    if column not in self.uniqueValues:
                        self.uniqueValues.append(column)
                # Get annotation name from header
                annotation = self.header[idx]
                self.data_subcollection[annotation]['values'].append(column)
            else:
                # If column isn't a annotation value, it's a cell name
                self.cell_names.append(column)

    def create_documents(self, file_path, file_id, study_accession):
        """Creats top level documents for Cell Metadata data structure"""
        documents = {}

        # Each annotation value has a top level document
        for value in self.header[1:]:
            # Copy document model so memory references are different
            copy_of_doc_model = copy.copy({
                'name': value,
                'study_accession': study_accession,
                'source_file_name': file_path.strip("."),
                'source_file_type': 'metadata',
                'annotation_type': ['group', 'numeric'],
                'file_id': file_id,
            })
            documents[value] = copy_of_doc_model
        return documents

    def create_subdocuments(self):
        """Creats subdocuments for each annotation """
        sub_documents = {}
        for value in self.header[1:]:
            # Copy subdocument model so memory references are different
            copy_of_subdoc_model = copy.copy({
                'cell_names': self.cell_names,
                'values': []
            })
            sub_documents[value] = copy_of_subdoc_model
        return sub_documents

    def get_collection_name(self):
        """Returns collection name"""
        return 'cell_metadata'

    def get_subcollection_name(self):
        """Returns sub-collection name"""
        return 'data'

    def chunk_subdocs(self):
