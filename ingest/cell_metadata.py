import argparse
import os
from typing import Dict, Generator, List, Tuple, Union

from ingest_pipeline_files import IngestFiles


class CellMetadata(IngestFiles):
    ALLOWED_FILE_TYPES = ['text/csv',
                          'text/plain', 'text/tab-separated-values']

    def __init__(self, file_path):
        IngestFiles.__init__(self, file_path, self.ALLOWED_FILE_TYPES, file_id, study_accession,
                             annotation_type, )
        self.header = self.file.readline().rstrip('\n').split(',')
        self.metadata_types = self.file.readline().rstrip('\n').split(',')
        # unique values for group-based annotations
        self.uniqueValues = []
        # subdocument data
        self.values = []
        self.cell_names = []

    def transform(self, line: List[str]) -> None:
        value = line.rstrip('\n').split(',')
        name_idx = value('NAME')
        for idx, x in enumerate(value):
             # determine whether or not value needs to be cast as a float or not
            if isinstance(x, int):
                self.values.append(round(float(x), 3))
            else:
                self.values.append(x)
                # determine if a new unique value needs to be stored in values array
                if self.metadata_types[idx] == 'group' and x not in self.uniqueValues:
                    self.uniqueValues.append(x)

    def create_document(self, name, annotation_type):
        return {'name': name,
                'study_accession': study_accession,
                'source_file_name': self.source_file_name
                'annotation_type': annotation_type,
                'unique_values': self.unique_values,
                'file_id', file_id,
                'data': {
                    'cell_names': [],
                    'values': []
                }
                }
