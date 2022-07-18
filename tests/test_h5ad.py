""" test_de.py
    integration test to verify that de process generates expected output
"""

import unittest
import sys
import hashlib
import os
import glob
from unittest.mock import patch
import scanpy as sc


sys.path.append("../ingest")
from ingest_files import IngestFiles
from h5ad import H5adIngestor


class TestH5adIngestor(unittest.TestCase):
    def test_known_good_obtain_adata(self):
        good_input = H5adIngestor(
            "../tests/data/h5ad/test.h5ad",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
        )
        self.assertTrue(
            good_input.validate(), "expect known good file to open with scanpy"
        )

    def test_known_bad_obtain_adata(self):
        bad_input = H5adIngestor(
            "../tests/data/h5ad/bad.h5",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
        )
        # passing obtain_data function to assertRaises using lambda
        # otherwise bad_input.obtain_data() is evaluated and triggers
        # an exception before assertRaises gets called
        self.assertRaises(ValueError, lambda: bad_input.obtain_adata())
        self.assertFalse(bad_input.validate())

