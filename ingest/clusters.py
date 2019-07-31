
import copy
import json
import re
import time
from typing import Dict, Generator, List, Tuple, Union

from ingest_pipeline_files import IngestFiles


class Clusters(IngestFiles):

    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']
    MAX_THRESHOLD = 100_000
    SUBSAMPLE_THRESHOLDS = [MAX_THRESHOLD, 20000, 10000, 1000]

    def __init__(self, file_path, *, name: str, study_accession: str = "", domain_ranges=None):

        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.header = self.split_line(self.file.readline().rstrip('\n'))
        # Second line in cluster is metadata_type
        self.metadata_types = self.split_line(
            self.file.readline().rstrip('\n'))
        self.uniqueValues = []
        self.source_file_name = file_path.strip(".")
        self.source_file_type = 'cluster'
        # self.points = amount of rows
        self.top_level_doc = {
            'name': name,
            'study_accession': study_accession,
            'source_file_name': file_path.strip("."),
            'source_file_type': 'cluster',
            'domain_ranges': domain_ranges,
            'points': self.amount_of_lines,
        }
        self.annotation_subdocs = self.create_annotation_subdocs()
        # sample size needs to be smaller than amount of points
        # for each cell annotations if there are any

    def transform(self, rows):
        """ Add data from cluster files into annotation subdocs in cluster data model"""
        data = self.split_line(rows)
        for idx, column in enumerate(data):
            if idx != 0:
                if self.metadata_types[idx].lower() == 'numeric':
                    column = round(float(column), 3)
                elif self.metadata_types[idx].lower() == 'group':
                    if column not in self.uniqueValues:
                        self.uniqueValues.append(column)
            annotation = self.header[idx].lower()
            # perform a shallow copy
            annotation_value = copy.copy(
                self.annotation_subdocs[annotation]['value'])
            annotation_value.append(column)
            self.annotation_subdocs[annotation]['value'] = annotation_value

    def create_annotation_subdocs(self):
        """Creates annotation_subdocs"""
        annotation_subdocs = dict.fromkeys(self.header)
        for value in self.header:
            value = value.lower()
            if value == 'name':
                annotation_subdocs[value] = self.create_metadata_subdoc(
                    'text', 'cells')
            elif value in ('x', 'y', 'z'):
                annotation_subdocs[value] = self.create_metadata_subdoc(
                    value, 'coordinates')
            else:
                annotation_subdocs[value] = self.create_metadata_subdoc(
                    value, 'annotations')
        return annotation_subdocs

    def get_collection_name(self):
        return 'clusters'

    def get_subcollection_name(self):
        return 'data'

    def create_metadata_subdoc(self, name, header_value_type, *, value=[], subsampled_annotation=None):
        return {
            'name': name,
            'array_index': 0,
            'value': value,
            'array_type': header_value_type,
            'subsampled_annotation': subsampled_annotation,
            'subsamp_threashold': "",
        }
