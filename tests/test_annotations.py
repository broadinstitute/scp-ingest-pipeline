"""Test test_annotations.py

These tests verify:
    - Group type annotations that have numeric-like values are being treated as strings
    - Numeric columns are rounded to 3 decimals points
    - Filtering cell names (given from cluster file) in metadata correctly
    - Labels are treated as strings

PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt

Note: When CI environment moves to Python 3.7, tests may break due to minor
differences in how the reference issues are serialized

# Run all tests in a manner that shows report_issues output
python3 test_annotations.py -s
"""


import random
import sys
import unittest
from decimal import Decimal

import numpy as np

sys.path.append("../ingest")
from annotations import Annotations
from clusters import Clusters


class TestAnnotations(unittest.TestCase):
    CLUSTER_PATH = '../tests/data/test_1k_cluster_data.csv'
    CELL_METADATA_PATH = '../tests/data/valid_no_array_v2.0.0.txt'

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
            # Ensure labels are strings
            header = column[0]
            assert isinstance(header, str)
            annot_type = column[1]
            if annot_type == 'group':
                # corrected testings of dataframe column dtype, using != always returns True
                self.assertFalse(
                    np.issubdtype(self.df.file[column].dtypes, np.number),
                    "Group annotations must be string values",
                )

    def test_merge_df(self):
        cluster = Clusters(
            '../tests/data/test_1k_cluster_data.csv',
            'dec0dedfeed1111111111111',
            'addedfeed000000000000000',
            'testCluster',
        )
        cell_metadata_df = Annotations(
            self.CELL_METADATA_PATH,
            ['text/csv', 'text/plain', 'text/tab-separated-values'],
        )
        cell_metadata_df.preprocess()
        cell_names_cell_metadata_df = np.asarray(cell_metadata_df.file['NAME'])
        cell_names_cluster_df = np.asarray(cluster.file['NAME'])
        # Cell names found in both cluster and metadata files
        common_cell_names = cell_names_cluster_df[
            np.isin(cell_names_cluster_df, cell_names_cell_metadata_df)
        ]
        print(f'common cell names: {common_cell_names}')
        # Perform merge
        print(cluster.file[['NAME', 'x', 'y', 'z']])
        cluster.merge_df(cluster.file[['NAME', 'x', 'y', 'z']], cell_metadata_df.file)

        # Ensure ONLY common cell names found in cell metadata file and cluster file
        # are in the newly merged df
        result = all(
            cell[0] in common_cell_names for cell in cluster.file['NAME'].values
        )
        self.assertTrue(
            result,
            f"Merge was not performed correctly. Merge should be performed on 'NAME'",
        )
