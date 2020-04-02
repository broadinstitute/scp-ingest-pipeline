"""Test test_annotations.py

These tests verify:
    - Group type annotations that have numeric-like values are being treated as strings
    - Numeric columns are rounded to 3 decimals points
    - Filtering cell names (given from cluster file) in metadata correctly

PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt

Note: When CI environment moves to Python 3.7, tests may break due to minor
differences in how the reference issues are serialized

# Run all tests in a manner that shows report_issues output
python3 test_ingest_files.py -s
"""


import sys
import unittest
from decimal import Decimal
import numpy as np
import random

sys.path.append("../ingest")

from annotations import Annotations


class TestAnnotations(unittest.TestCase):
    CLUSTER_PATH = '../tests/data/test_1k_cluster_data.csv'
    CELL_METADATA_PATH = '../tests/data/valid_no_array_v2.0.0.tsv'

    EXPONENT = -3

    def setUp(self):
        self.df = Annotations(
            self.CLUSTER_PATH, ['text/csv', 'text/plain', 'text/tab-separated-values']
        )

    def test_round(self):
        # Pick a random number between 1 and amount of lines in file
        ran_num = random.randint(1, 2000)
        self.df.preprocess()
        for column in self.df.file.columns:
            annot_type = column[1]
            if annot_type == 'numeric':
                value = str(self.df.file[column][ran_num])
                print(Decimal(value).as_tuple().exponent)
                assert (
                    abs(Decimal(value).as_tuple().exponent) >= self.EXPONENT
                ), "Numbers did not round to 3 or less decimals places"

    def test_group_annotations(self):
        self.df.preprocess()
        for column in self.df.file.columns:
            annot_type = column[1]
            if annot_type == 'group':
                assert (
                    self.df.file[column].dtype != np.number
                ), "Group annotations must be string values"

    def test_merge_df(self):
        self.df.preprocess()
        cell_metadata_df = Annotations(
            self.CELL_METADATA_PATH,
            ['text/csv', 'text/plain', 'text/tab-separated-values'],
        )
        cell_metadata_df.preprocess()
        cell_names_cell_metadata_df = np.asarray(cell_metadata_df.file['NAME'])
        cell_names_cluster_df = np.asarray(self.df.file['NAME'])

        # Cell names found in both cluster and metadata files
        common_cell_names = cell_names_cluster_df[
            np.isin(cell_names_cluster_df, cell_names_cell_metadata_df)
        ]
        # Perform merge
        self.df.merge_df(self.df.file[['NAME', 'X', 'Y', 'Z']], cell_metadata_df.file)

        # Ensure ONLY common cell names found in cell metadata file and cluster file
        # are in the newly merged df
        result = all(
            cell[0] in common_cell_names for cell in self.df.file['NAME'].values
        )
        self.assertTrue(
            result,
            f"Merge was not performed correctly. Merge should be performed on 'NAME'",
        )
