""" test_anndata.py
    verify basic AnnData validation works as expected
"""

import unittest
import sys

sys.path.append("../ingest")
from anndata_ import AnnDataIngestor


class TestAnnDataIngestor(unittest.TestCase):
    def test_minimal_valid_anndata(self):
        good_input = AnnDataIngestor(
            "../tests/data/anndata/test.h5ad",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
        )
        self.assertTrue(
            good_input.validate(), "expect known good file to open with scanpy"
        )

    def test_truncated_anndata(self):
        truncated_input = AnnDataIngestor(
            "../tests/data/anndata/bad.h5",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
        )
        # passing obtain_data function to assertRaises using lambda
        # otherwise truncated_input.obtain_data() is evaluated and triggers
        # an exception before assertRaises gets called
        self.assertRaisesRegex(
            ValueError,
            "Scanpy cannot read file, \"../tests/data/anndata/bad.h5\".",
            lambda: truncated_input.obtain_adata(),
        )
        self.assertFalse(truncated_input.validate())

    def test_input_bad_suffix(self):
        bad_input = AnnDataIngestor(
            "../tests/data/anndata/bad.foo",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
        )
        # passing obtain_data function to assertRaises using lambda
        # otherwise bad_input.obtain_data() is evaluated and triggers
        # an exception before assertRaises gets called
        self.assertRaisesRegex(
            ValueError,
            "File type not detected for ../tests/data/anndata/bad.foo, expected file endings are: .h5ad .h5 .hdf5",
            lambda: bad_input.obtain_adata(),
        )
        self.assertFalse(bad_input.validate())

