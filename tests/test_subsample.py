"""Tests subsampling

These tests verify:
    - Binning correctly
    - Subsampling correctly

PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt

Note: When CI environment moves to Python 3.7, tests may break due to minor
differences in how the reference issues are serialized

# Run all tests in a manner that shows report_issues output
python3 test_subsample.py -s

"""

import sys
import unittest
import numpy as np

sys.path.append("../ingest")

from subsample import SubSample


class TestSubsample(unittest.TestCase):
    AMOUNT_OF_NUMERIC_BINS = 20
    CLUSTER_PATH = '../tests/data/test_1k_cluster_data.csv'

    def setUp(self):
        self.subsample_obj = SubSample(self.CLUSTER_PATH)

    def test_subsample(self):
        for data in self.subsample_obj.subsample():
            header_value = data[1]
            annot_name = data[1][0].lower()
            annot_type = data[1][1]
            subsampled_data = data[0]
            if annot_type == 'group':
                amount_of_bins = len(self.subsample_obj.file[header_value].unique())
            else:
                amount_of_bins = self.AMOUNT_OF_NUMERIC_BINS
            for key_value in subsampled_data.items():
                # values from current column in orginal df
                print(f"Current annotation is {key_value[0]}")
                orig_vals_lists = np.array_split(
                    np.concatenate(self.subsample_obj.file[key_value[0]].values),
                    amount_of_bins,
                )

                subsampled_vals_lists = np.array_split(key_value[1], amount_of_bins)

                for (orig_val_list, subsampled_val_list) in zip(
                    orig_vals_lists, subsampled_vals_lists
                ):

                    result = np.isin(subsampled_val_list, orig_val_list)

                    bad_values = np.isin(
                        subsampled_val_list, orig_val_list, invert=True
                    )

                    self.assertTrue(
                        all(result),
                        f"Incorrect subsampling. Values that should not have been in subsampled for {key_value[0]} in {annot_name}:  {subsampled_val_list[bad_values]}",
                    )

    def test_bin(self):
        for bin_data in map(self.subsample_obj.bin, self.subsample_obj.columns):
            bins = bin_data[0]
            column_name = bin_data[1]
            annot_type = bin_data[1][1]
            if annot_type == 'group':
                expected_unique_groups = self.subsample_obj.file[column_name].unique()

                print(f"Expected values {expected_unique_groups}")
                # Current unique values that has a bin
                current_unique_groups = list(bins.keys())
                print(current_unique_groups)
                result = all(
                    group in expected_unique_groups for group in current_unique_groups
                )
                # Group annotations should have a bin for each unique value
                self.assertTrue(
                    result,
                    f"Each unique value should have its own bin for a given column. Expected bins for {self.subsample_obj.file[column_name].unique()}",
                )
            else:
                # Numeric annotations should have 20 bins
                self.assertEqual(
                    len(bins.keys()),
                    self.AMOUNT_OF_NUMERIC_BINS,
                    "Metadata validation issues do not match reference issues",
                )
