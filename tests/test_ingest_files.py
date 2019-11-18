"""Test ingest_files.py
These tests verify:
    - Error is thrown for a missing file
PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
Note: When CI environment moves to Python 3.7, tests may break due to minor
differences in how the reference issues are serialized

# Run all tests in a manner that shows report_issues output
python3 test_ingest_files.py -s
"""


import sys
import unittest

sys.path.append("../ingest")
from ingest_files import IngestFiles


class TestIngestFiles(unittest.TestCase):
    def test_ingest_missing_file(self):
        """Should throw error for missing local file
        """
        with self.assertRaises(OSError):
            IngestFiles(
                '/this/file/does/not_exist.txt',
                ['text/csv', 'text/plain', 'text/tab-separated-values'],
                open_as='dataframe',
            )
