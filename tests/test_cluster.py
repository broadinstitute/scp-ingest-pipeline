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
