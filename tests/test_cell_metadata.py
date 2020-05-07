import sys
import unittest

sys.path.append("../ingest")
from cell_metadata import CellMetadata

class TestCellMetadata(unittest.TestCase):

    def test_validate_header_for_coordinate_values_false(self):
        """Ensures validate_header_for_coordinate_values returns false when
        coordintate value is in metadata file
         """
        cm = CellMetadata(
            '../tests/data/metadata_bad_contains_coordinates.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            study_accession='SCP2',
            tracer=None,
        )
        self.assertFalse(cm.validate_header_for_coordinate_values())

    def test_validate_header_for_coordinate_values_true(self):
        """Ensures validate_header_for_coordinate_values returns true when
        coordintate value is not in metadata file
         """
        cm = CellMetadata(
            '../tests/data/metadata_example.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            study_accession='SCP2',
            tracer=None,
        )
        self.assertTrue(cm.validate_header_for_coordinate_values())
