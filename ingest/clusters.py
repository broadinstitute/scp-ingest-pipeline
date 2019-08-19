import copy
import unicodedata
from typing import Dict, Generator, List, Tuple, Union

from ingest_files import IngestFiles


class Clusters(IngestFiles):

    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']
    COLLECTION_NAME = 'clusters'
    SUBCOLLECTION_NAME = 'data'

    def __init__(self, file_path, *, name: str = None, study_accession: str = "", domain_ranges=None):

        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.header = self.get_next_line(increase_line_count=False)
        # Second line in cluster is metadata_type
        self.metadata_types = self.get_next_line(increase_line_count=False)
        self.unique_values = dict.fromkeys(self.header[1:], [])
        self.source_file_type = 'cluster'
        self.has_z = 'z' in self.header
        self.top_level_doc = {
            'name': name,
            'cluster_type': '3d' if self.has_z else '2d',
            'cell_annotations': [],
            'study_accession': study_accession,
            'source_file_type': 'cluster',
            'domain_ranges': domain_ranges,
            'points': self.amount_of_lines,
        }
        self.cluster_subdocs = self.return_cluster_subdocs()

    def update_points(self):
        self.top_level_doc['points'] = self.amount_of_lines

    def transform(self, row):
        """ Add data from cluster files into annotation subdocs in cluster data model"""

        for idx, column in enumerate(row):
            annotation = self.header[idx].casefold()

            # first index is cell name don't need to check annot type
            if idx != 0:
                if self.metadata_types[idx].casefold() == 'numeric':
                    column = round(float(column), 3)
                elif self.metadata_types[idx].casefold() == 'group':
                    if column not in self.unique_values[annotation]:
                        self.unique_values[annotation].append(column)
            # perform a shallow copy
            annotation_value = copy.copy(
                self.cluster_subdocs[annotation]['value'])
            annotation_value.append(column)
            self.cluster_subdocs[annotation]['value'] = annotation_value

    def return_cluster_subdocs(self):
        """Creates cluster_subdocs"""
        cluster_subdocs = {}
        for annot_name in self.header:
            value = annot_name.casefold()
            if value == 'name':
                cluster_subdocs[value] = self.create_cluster_subdoc(
                    'text', 'cells')
            elif value in ('x', 'y', 'z'):
                cluster_subdocs[value] = self.create_cluster_subdoc(
                    value, 'coordinates')
            else:
                cluster_subdocs[value] = self.create_cluster_subdoc(
                    annot_name, 'annotations')
        return cluster_subdocs

    def create_cluster_subdoc(self, annot_name, header_value_type, *, value=[], subsample_annotation=None):
        """Returns cluster subdoc"""
        return {
            'name': annot_name,
            'array_index': 0,
            'value': value,
            'array_type': header_value_type,
            'subsample_annotation': subsample_annotation,
            'subsample_threshold': '',
        }

    def can_subsample(self):
        # TODO: Add more subsample validations
        if self.has_z:
            return len(self.header) > 4
        else:
            return len(self.header) > 3
