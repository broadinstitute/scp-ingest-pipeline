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
    CLUSTER_PATH = "../tests/data/test_1k_cluster_data.csv"
    CELL_METADATA_PATH = "../tests/data/valid_no_array_v2.0.0.txt"

    EXPONENT = -3

    def setUp(self):
        self.df = Annotations(
            self.CLUSTER_PATH, ["text/csv", "text/plain", "text/tab-separated-values"]
        )

    def test_create_columns(self):
        header = ["Intensity", "donor_id", "species__ontology_label"]
        annotatiion_types = ["numeric", "group", "group"]
        colums = Annotations.create_columns(header, annotatiion_types)
        expected = [
            ("Intensity", "numeric"),
            ("donor_id", "group"),
            ("species__ontology_label", "group"),
        ]
        self.assertEqual(colums, expected)

    def test_duplicate_headers(self):
        """Annotation headers should not contain duplicate values
        """
        dup_headers = Annotations(
            "../tests/data/dup_headers_v2.0.0.tsv",
            ["text/csv", "text/plain", "text/tab-separated-values"],
        )

        self.assertFalse(
            dup_headers.validate_unique_header(),
            "Duplicate headers should fail format validation",
        )

        with self.assertRaises(ValueError):
            dup_headers.preprocess()

    def test_set_dtypes(self):
        headers = ["NAME", "", ""]
        annot_types = ["TYPE", "GROUP", "numeric"]

    def test_header_format(self):
        """Header rows of metadata file should conform to standard
        """
        error_headers = Annotations(
            "../tests/data/error_headers_v2.0.0.tsv",
            ["text/csv", "text/plain", "text/tab-separated-values"],
        )

        self.assertFalse(
            error_headers.validate_header_keyword(),
            "Missing NAME keyword should fail format validation",
        )

        self.assertFalse(
            error_headers.validate_type_keyword(),
            "Missing TYPE keyword should fail format validation",
        )

        self.assertFalse(
            error_headers.validate_type_annotations(),
            "Invalid type annotations should fail format validation",
        )

    def test_low_mem_artifact(self):
        # pandas default of low_memory=True allows internal chunking during parsing
        # causing inconsistent dtype coercion artifact for larger annotation files

        lmtest = Annotations(
            "../tests/data/low_mem_unit.txt",
            ["text/csv", "text/plain", "text/tab-separated-values"],
        )
        lmtest.preprocess()

        # when low memory=True, the first row in the file would be in the first chunk
        # and the numeric value was not properly coerced to become a string
        assert isinstance(
            lmtest.file["mixed_data"]["group"][0], str
        ), "numeric value should be coerced to string"

        # Per SCP-2545 NA values become strings for group annotations.
        print(lmtest.file["mixed_data"]["group"][2])
        print(type(lmtest.file["mixed_data"]["group"][2]))
        assert isinstance(
            lmtest.file["mixed_data"]["group"][2], str
        ), "expect empty cell conversion to NaN is string for group annotation"

        # numeric value in second chunk should still properly be coerced to string type
        assert isinstance(
            lmtest.file["mixed_data"]["group"][32800], str
        ), "numeric value should be coerced to string"

    def test_round_numeric_annotations(self):
        # Pick a random number between 1 and amount of lines in file
        ran_num = random.randint(1, 2000)
        self.df.preprocess()
        for column in self.df.file.columns:
            annot_type = column[1]
            if annot_type == "numeric":
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
            if annot_type == "group":
                # corrected testings of dataframe column dtype, using != always returns True
                self.assertFalse(
                    np.issubdtype(self.df.file[column].dtypes, np.number),
                    "Group annotations must be string values",
                )

    def test_merge_df(self):
        cluster = Clusters(
            "../tests/data/test_1k_cluster_data.csv",
            "dec0dedfeed1111111111111",
            "addedfeed000000000000000",
            "testCluster",
        )
        cell_metadata_df = Annotations(
            self.CELL_METADATA_PATH,
            ["text/csv", "text/plain", "text/tab-separated-values"],
        )
        cell_metadata_df.preprocess()
        cell_names_cell_metadata_df = np.asarray(cell_metadata_df.file["NAME"])
        cell_names_cluster_df = np.asarray(cluster.file["NAME"])
        # Cell names found in both cluster and metadata files
        common_cell_names = cell_names_cluster_df[
            np.isin(cell_names_cluster_df, cell_names_cell_metadata_df)
        ]
        print(f"common cell names: {common_cell_names}")
        # Perform merge
        print(cluster.file[["NAME", "x", "y", "z"]])
        cluster.merge_df(cluster.file[["NAME", "x", "y", "z"]], cell_metadata_df.file)

        # Ensure ONLY common cell names found in cell metadata file and cluster file
        # are in the newly merged df
        result = all(
            cell[0] in common_cell_names for cell in cluster.file["NAME"].values
        )
        self.assertTrue(
            result,
            f"Merge was not performed correctly. Merge should be performed on 'NAME'",
        )

    def test_get_cell_names(self):
        import pandas as pd

        expected_cell_names = [
            "CELL_0001",
            "  CELL_0002",
            "  CELL_0003",
            "  CELL_0004",
            "  CELL_0005",
            "  CELL_0006",
            "  CELL_0007",
            "  CELL_0008",
            "  CELL_0009",
            " CELL_00010",
            " CELL_00011",
            " CELL_00012",
            " CELL_00013",
            " CELL_00014",
            " CELL_00015",
            " CELL_00016",
            " CELL_00017",
            " CELL_00018",
            " CELL_00019",
            " CELL_00020",
        ]
        column_names = [
            ("NAME", "TYPE"),
            ("Cluster", "group"),
            ("Sub-Cluster", "group"),
            ("Average Intensity", "numeric"),
        ]
        index = pd.MultiIndex.from_tuples(column_names)

        df = pd.read_csv(
            "../tests/data/metadata_example.txt", sep="\t", names=index, skiprows=2
        )
        cells = Annotations.get_cell_names(df)
        self.assertEqual(cells, expected_cell_names)
