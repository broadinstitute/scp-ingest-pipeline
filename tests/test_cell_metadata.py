import sys
import unittest
from unittest.mock import MagicMock, patch

sys.path.append("../ingest")
from cell_metadata import CellMetadata
from annotations import Annotations


class TestCellMetadata(unittest.TestCase):
    def setUp(self):
        # This cluster file should be sampled at 1k
        self.cm = CellMetadata(
            '../tests/data/metadata_example.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            study_accession= 'SCP2',
            tracer=None,
        )

    def test_validate_gene_keyword_false(self):
        "Validates validate_gene_keyword() returns false correctly"
        cm = CellMetadata(
            '../tests/data/metadata_bad.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            study_accession= 'SCP2',
            tracer=None,
        )
        self.assertFalse(cm.validate_header_for_coordinate_values())

    def test_validate_gene_keyword_true(self):
        "Validates validate_gene_keyword() returns false correctly"
        self.assertTrue(self.cm.validate_header_for_coordinate_values())

    @patch("cell_metadata.validate_header_for_coordinate_values", MagicMock(return_value="True"))
    @patch("annotations.validate_format", side_effect=False)
    def test_is_valid_format_false(self, mock_validate_format, mock_validate_header_for_coordinate_values):
        "Checks to see if is_valid_format returns false"
        assertFalse(self.cm.is_valid_format())

    @patch("CellMetadata.validate_header_for_coordinate_values", MagicMock(return_value="True"))
    @patch("Annotations.validate_format", side_effect=True)
    def test_is_valid_format_true(self):
        "Validates validate_gene_keyword() returns true correctly"
        self.assertTrue(self.cm.is_valid_format())
