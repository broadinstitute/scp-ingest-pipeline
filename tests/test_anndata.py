""" test_anndata.py
    verify basic AnnData validation works as expected
"""

import unittest
import sys
import os
from unittest.mock import patch

sys.path.append("../ingest")
from anndata_ import AnnDataIngestor
from ingest_files import IngestFiles


class TestAnnDataIngestor(unittest.TestCase):
    @staticmethod
    def setup_class(self):
        filepath_valid = "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad"
        filepath_invalid = "../tests/data/anndata/bad.h5"
        self.study_id = "addedfeed000000000000000"
        self.study_file_id = "dec0dedfeed0000000000000"
        self.valid_args = [filepath_valid, self.study_id, self.study_file_id]
        self.invalid_args = [filepath_invalid, self.study_id, self.study_file_id]
        self.cluster_name = 'X_tsne'
        self.valid_kwargs = {'obsm_keys': [self.cluster_name]}
        self.anndata_ingest = AnnDataIngestor(*self.valid_args, **self.valid_kwargs)
        self.output_filename = (
            f"{self.study_file_id}.{self.cluster_name}.cluster.anndata_segment.tsv"
        )

    def teardown_method(self, _):
        if os.path.isfile(self.output_filename):
            os.remove(self.output_filename)

    def test_minimal_valid_anndata(self):
        self.assertTrue(
            self.anndata_ingest.validate(), "expect known good file to open with scanpy"
        )

    def test_truncated_anndata(self):
        truncated_input = AnnDataIngestor(*self.invalid_args)
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
            "../tests/data/anndata/bad.foo", self.study_id, self.study_file_id
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

    def test_set_output_filename(self):
        cluster_name = "dec0dedfeed0000000000000.X_Umap"
        self.assertEqual(
            AnnDataIngestor.set_clustering_filename(cluster_name),
            "dec0dedfeed0000000000000.X_Umap.cluster.anndata_segment.tsv",
        )

    def test_generate_cluster_header(self):
        self.anndata_ingest.generate_cluster_header(
            self.anndata_ingest.obtain_adata(), self.cluster_name, self.study_file_id
        )
        with open(self.output_filename) as header_file:
            header = header_file.readline().split("\t")
            self.assertEqual(
                ['NAME', 'X', "Y\n"], header, "did not find expected headers"
            )

    def test_generate_cluster_type_declaration(self):
        self.anndata_ingest.generate_cluster_type_declaration(
            self.anndata_ingest.obtain_adata(), self.cluster_name, self.study_file_id
        )
        with open(self.output_filename) as header_file:
            header = header_file.readline().split("\t")
            self.assertEqual(
                ['TYPE', 'numeric', "numeric\n"],
                header,
                "did not find expected headers",
            )

    def test_generate_cluster_body(self):
        self.anndata_ingest.generate_cluster_body(
            self.anndata_ingest.obtain_adata(), self.cluster_name, self.study_file_id
        )
        with open(self.output_filename) as cluster_body:
            line = cluster_body.readline().split("\t")
            expected_line = ['AAACATACAACCAC-1', '16.009954', "-21.073845\n"]
            self.assertEqual(
                expected_line,
                line,
                'did not get expected coordinates from cluster body',
            )

    def test_get_files_to_delocalize(self):
        files = AnnDataIngestor.clusterings_to_delocalize(
            self.valid_kwargs, self.study_file_id
        )
        expected_files = [self.output_filename]
        self.assertEqual(expected_files, files)

    def test_delocalize_files(self):
        # just create header, no reason to run full extract
        self.anndata_ingest.generate_cluster_header(
            self.anndata_ingest.obtain_adata(), self.cluster_name, self.study_file_id
        )
        with patch('ingest_files.IngestFiles.delocalize_file'):
            AnnDataIngestor.delocalize_file(
                "gs://fake_bucket",
                self.study_id,
                AnnDataIngestor.clusterings_to_delocalize(
                    self.valid_kwargs, self.study_file_id
                ),
            )
            self.assertEqual(
                IngestFiles.delocalize_file.call_count,
                1,
                "expected 1 call to delocalize output files",
            )
