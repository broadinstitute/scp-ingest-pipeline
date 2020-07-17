import sys
import unittest

sys.path.append("../ingest")
from cell_metadata import CellMetadata


class TestCellMetadata(unittest.TestCase):
    def test_validate_header_for_coordinate_values_false(self):
        """Ensures validate_header_for_coordinate_values returns false when
        coordinate value is in metadata file
        Note: cluster has similar set of tests
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
        coordinate value is not in metadata file
         """
        cm = CellMetadata(
            '../tests/data/metadata_example.txt',
            '5d276a50421aa9117c982845',
            '5dd5ae25421aa910a723a337',
            study_accession='SCP2',
            tracer=None,
        )
        self.assertTrue(cm.validate_header_for_coordinate_values())

    def test_metadata_type_inference(self):
        """Confirm consistency of type inference behavior
        in instantiated data frame
        """
        cm = CellMetadata(
            '../tests/data/metadata_NA.txt',
            'addedfeed000000000000000',
            'dec0dedfeed1111111111111',
            study_accession='SCPtest',
            tracer=None,
        )

        # integers, empty cell and string as inputs for numeric annotation
        assert isinstance(
            cm.file['NA_i_n_s__grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        # integers and empty cell as inputs for numeric annotation
        assert isinstance(
            cm.file['NA_i_n_grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        # floats, empty cell and string as inputs for numeric annotation
        assert isinstance(
            cm.file['NA_f_n_s__grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        # floats and empty cell as inputs for numeric annotation
        assert isinstance(
            cm.file['NA_f_n_grp']['group'][3], str
        ), "empty cell -> NaN, expect coercion to string"

        # integers, empty cell and string as inputs for group annotation
        assert isinstance(
            cm.file['NA_i_n_s__num']['numeric'][3], float
        ), "empty cell -> NaN that remains float (not coerced)"

        # floats, empty cell and string as inputs for group annotation
        assert isinstance(
            cm.file['NA_f_n_s__num']['numeric'][3], float
        ), "empty cell -> NaN that remains float (not coerced)"

        self.assertFalse(
            cm.validate_numeric_annots(),
            'numeric annotations supplied with strings should be invalid',
        )
