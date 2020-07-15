import sys
import unittest

sys.path.append("../ingest")
from clusters import Clusters


class TestCellMetadata(unittest.TestCase):
    def test_validate_header_for_coordinate_values_false(self):
        """Ensures validate_header_for_coordinate_values returns false when
         coordinate is missing in header
        """
        cluster = Clusters(
            '../tests/data/cluster_bad_missing_coordinate.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            'testCluster',
        )
        self.assertFalse(cluster.validate_header_for_coordinate_values())

    def test_validate_header_for_coordinate_values_true(self):
        """Ensures validate_header_for_coordinate_values returns true when
        coordintate value is in cluster file
         """
        cluster = Clusters(
            '../tests/data/cluster_example.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            'testCluster',
        )
        self.assertTrue(cluster.validate_header_for_coordinate_values())

    def test_non_coordinate_case_integrity(self):
        # if headers of this file change, the reference for the assert needs to be updated
        cluster = Clusters(
            '../tests/data/test_1k_cluster_data.csv',
            'dec0dedfeed1111111111111',
            'addedfeed000000000000000',
            'testCluster',
        )
        assert cluster.headers == [
            'NAME',
            'x',
            'y',
            'z',
            'CLUSTER',
            'SUBCLUSTER',
            '1',
        ], 'cluster instantiation should only downcase coordinate header columns'

    def test_cluster_type_inference(self):
        """Confirm consistency of type inference behavior
        in instantiated data frame
        Note: metadata has similar set of tests
        """
        cluster = Clusters(
            '../tests/data/cluster_NA.txt',
            'addedfeed000000000000000',
            'dec0dedfeed1111111111111',
            'testCluster',
        )

        assert isinstance(
            cluster.file['NA_i_n_s__grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        assert isinstance(
            cluster.file['NA_f_n_s__grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        assert isinstance(
            cluster.file['NA_i_n_s__num']['numeric'][3], float
        ), "empty cell -> NaN that remains float (not coerced)"

        assert isinstance(
            cluster.file['NA_f_n_s__num']['numeric'][3], float
        ), "empty cell -> NaN that remains float (not coerced)"

        self.assertFalse(
            cluster.validate_numeric_annots(),
            'numeric annotations supplied with strings should be invalid',
        )
