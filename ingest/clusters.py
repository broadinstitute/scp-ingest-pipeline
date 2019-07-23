
import time
from typing import Dict, Generator, List, Tuple, Union

from ingest_pipeline_files import IngestFiles


class Clusters(IngestFiles):

    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']
    MAX_THRESHOLD = 100_000
    SUBSAMPLE_THRESHOLDS = [MAX_THRESHOLD, 20000, 10000, 1000]

    def __init__(self, file_path, *, name: str = "", study_accession: str = "", points=None, domain_ranges=None):
        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES)
        self.header = self.file.readline().rstrip('\n').split(',')
        self.metadata_types = self.file.readline().rstrip('\n').split(',')
        self.uniqueValues = []
        self.source_file_name = file_path.strip(".")
        self.source_file_type = 'cluster'
        # self.points = amount of rows
        self.top_level_doc = {
        'name': name,
        'study_accession': study_accession,
        'source_file_name': source_file_name,
        'source_file_type': source_file_type,
        'domain_ranges': domain_ranges,
        'points': points,
        }
        self.annotation_subdocs = self.create_annotation_subdocs()
        # sample size needs to be smaller than amount of points
        # for each cell annotations if there are any

    def is_3d(self):
        return 'z' in self.header

    def transform_cluster_file(self, rows):
        data = rows.rstrip('\n').split(',')
        amount_of_rows = None
        for idx, row_value in enumerate(data):
            if metadata_types[idx] == 'numeric':
                round(float(row_value), 3)
            self.annotation_subdocs[idx].value.append(row_value)
            if metadata_types[idx] == 'group':
                if row_value not in self.uniqueValues:
                    self.uniqueValues.append(row_value)
            amount_of_rows = row_value
        # create subsampled data_arrays for visualization
        for amount_of_rows > sample in self.SUBSAMPLE_THRESHOLDS:

    def create_sub_sampled(self):

    def create_annotation_subdocs(self):
        types = iter(self.metadata_types)
        annotation_subdocs = {}
        next(types)
        for idx, value in enumerate(types):
            if value == 'name':
                annotation_subdocs.update(
                    value, self.create_metadata_subdoc('text', idx, 'cells'))
            elif value in ('x', 'y', 'z'):
                annotation_subdocs.update(
                    value, create_metadata_subdoc(value, idx, 'coordinates'))
            else:
                annotation_subdocs.update(
                    value, create_metadata_subdoc(value, idx, 'annotations'))
        return annotation_subdocs

    def create_metadata_subdoc(self, name, array_idx, array_type, *, value=[]], subsampled_annotation=None):
        return {
            'name': name,
            'array_index': array_idx,
            'value': value,
            'array_type': "",
            'subsampled_annotation': subsampled_annotation,
            'subsamp_threashold': "",
        }
