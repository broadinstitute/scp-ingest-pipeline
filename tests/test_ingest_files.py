import unittest
import sys

sys.path.append("../ingest")
from ingest_files import IngestFiles


class TestIngestFiles(unittest.TestCase):
    def test_is_gzipped(self):
        zipped_path = "../tests/data/mtx/AB_toy_data_toy.unsorted_mtx.mtx.gz"
        self.assertTrue(IngestFiles.is_gzipped(zipped_path))
        not_zipped = "../tests/data/expression_matrix_bad_duplicate_gene.txt"
        self.assertFalse(IngestFiles.is_gzipped(not_zipped))
