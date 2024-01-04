import unittest
import sys
import string
import random
from unittest.mock import patch
from google.cloud import storage

sys.path.append("../ingest")
from ingest_files import IngestFiles


class TestIngestFiles(unittest.TestCase):
    # Test exception

    def raise_error(*args):
        raise Exception

    def test_is_gzipped(self):
        zipped_path = "../tests/data/mtx/unsorted_matrix.mtx.gz"

        self.assertTrue(IngestFiles.is_gzipped(zipped_path))
        gzip_extension = "../tests/data/mtx/unsorted_mtx.txt.gzip"
        self.assertTrue(IngestFiles.is_gzipped(gzip_extension))
        not_zipped = "../tests/data/expression_matrix_bad_duplicate_gene.txt"
        self.assertFalse(IngestFiles.is_gzipped(not_zipped))

        # Test exception
        with patch("gzip.open", side_effect=TestIngestFiles.raise_error):
            self.assertRaises(ValueError, IngestFiles.is_gzipped, zipped_path)

    def test_delocalize_file(self):
        """Tests writing a small file to a GCP bucket
        Checks the test location doesn't already have the test file
        writes the file, checks it now exists, then deletes the test file
        Due to auth issues, this test does not run in CI
        """
        dev_reference_bucket = "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982"
        file_to_delocalize = (
            "../tests/data/differential_expression/sparse/sparsemini_features.tsv"
        )
        path_segments = file_to_delocalize.split("/")
        file = path_segments[-1]
        alphanum = string.ascii_letters + string.digits
        random_string = ''.join(random.choice(alphanum) for _ in range(10))
        dest_name = f"_pytest/{file}_{random_string}"
        print(f"dest_name")
        storage_client = storage.Client()
        bucket_name = dev_reference_bucket[5:]
        bucket = storage_client.bucket(bucket_name)

        try:
            IngestFiles.delocalize_file(
                dev_reference_bucket, file_to_delocalize, dest_name
            )
            self.assertTrue(
                storage.Blob(bucket=bucket, name=dest_name).exists(storage_client)
            )
        finally:
            storage.Blob(bucket=bucket, name=dest_name).delete()
