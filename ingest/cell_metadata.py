import argparse
import os
from typing import Dict, Generator, List, Tuple, Union


class Cell_Metadata:
    def __init__(self, file_path):
        self.file = open(file_path, 'r')
        self.header = self.file.readline().rstrip('\n').split(',')
        self.metadata_types = self.file.readline().rstrip('\n').split(',')
        self.annotation_type = ['group', 'numeric']
        # unique values for group-based annotations
        self.uniqueValues = []
        # subdocument data
        self.values = []
        self.cell_names

    def extract(self):
        for line in self.file.readlines():
            yield line

    def transform(self, line: List[str]) - >None:
        value = line.rstrip('\n').split(','):
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
