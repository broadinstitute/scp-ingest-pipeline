"""Unit Test for cell metadata files"""
import io
import sys
import unittest

sys.path.append('../ingest')
from cell_metadata import CellMetadata


class TestCellMetadata(unittest.TestCase):
    def test_chunk_subdocuments(self):
        self.cell_metadata = CellMetadata(
            '../tests/data/test_chunking_cell_metadata.csv', 'SCP1', '1234abc'
        )
        self.cell_metadata.preproccess()
        # Captures print statements from chunk_subdocuments
        captured_output = io.StringIO()
        sys.stdout = captured_output
        for cellmetadatadatamodel in self.cell_metadata.transform():
            for subdoc in self.cell_metadata.chunk_subdocuments(
                'doc_name', 'fake/path/cell_metadata', cellmetadatadatamodel
            ):
                pass
        sys.stdout = sys.__stdout__
        # String represents outputs from chunk_subdocuments method
        #  number of bytes, end_index, starting index, number of bytes, end_index, starting index
        expected = '1048562, 17740, 0, 17739\n347983, 23593, 17740, 23593\n1048567, 16130, 0, 16129\n485212, 23593, 16130, 23593\n'

        self.assertEqual(expected, captured_output.getvalue())


if __name__ == "__main__":
    unittest.main()
