""" test_anndata.py
    verify basic AnnData validation works as expected
"""

import unittest
import sys
import glob
import os
import gzip
import time
from unittest.mock import patch

from test_expression_files import mock_expression_load

sys.path.append("../ingest")
from anndata_ import AnnDataIngestor
from ingest_files import IngestFiles
from expression_files.expression_files import GeneExpression


class TestAnnDataIngestor(unittest.TestCase):
    GeneExpression.load = mock_expression_load

    @staticmethod
    def setup_class(self):
        filepath_valid = "../tests/data/anndata/trimmed_compliant_pbmc3K.h5ad"
        filepath_invalid = "../tests/data/anndata/bad.h5"
        filepath_dup_feature = "../tests/data/anndata/dup_feature.h5ad"
        filepath_dup_cell = "../tests/data/anndata/dup_cell.h5ad"
        filepath_nan = "../tests/data/anndata/nan_value.h5ad"
        filepath_synthetic = "../tests/data/anndata/anndata_test.h5ad"
        self.study_id = "addedfeed000000000000000"
        self.study_file_id = "dec0dedfeed0000000000000"
        self.valid_args = [filepath_valid, self.study_id, self.study_file_id]
        self.invalid_args = [filepath_invalid, self.study_id, self.study_file_id]
        self.dup_feature_args = [
            filepath_dup_feature,
            self.study_id,
            self.study_file_id,
        ]
        self.dup_cell_args = [filepath_dup_cell, self.study_id, self.study_file_id]
        self.nan_value_args = [filepath_nan, self.study_id, self.study_file_id]
        self.synthetic_args = [filepath_synthetic, self.study_id, self.study_file_id]
        self.cluster_name = 'X_tsne'
        self.valid_kwargs = {'obsm_keys': [self.cluster_name]}
        self.anndata_ingest = AnnDataIngestor(*self.valid_args, **self.valid_kwargs)
        self.cluster_filename = f"h5ad_frag.cluster.{self.cluster_name}.tsv"
        self.metadata_filename = "h5ad_frag.metadata.tsv"

    def teardown_method(self, _):
        if os.path.isfile(self.cluster_filename):
            os.remove(self.cluster_filename)
        for f in glob.glob("h5ad_frag*.gz"):
            os.remove(f)

    def test_minimal_valid_anndata(self):
        self.assertTrue(
            self.anndata_ingest.basic_validation(),
            "expect known good file to open with scanpy",
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
        self.assertFalse(truncated_input.basic_validation())

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
        self.assertFalse(bad_input.basic_validation())

    def test_set_output_filename(self):
        cluster_name = "X_Umap"
        self.assertEqual(
            AnnDataIngestor.set_clustering_filename(cluster_name),
            "h5ad_frag.cluster.X_Umap.tsv",
        )

    def test_generate_cluster_header(self):
        self.anndata_ingest.generate_cluster_header(
            self.anndata_ingest.obtain_adata(), self.cluster_name
        )
        with open(self.cluster_filename) as header_file:
            header = header_file.readline().split("\t")
            self.assertEqual(
                ['NAME', 'X', "Y\n"], header, "did not find expected headers"
            )

    def test_generate_cluster_type_declaration(self):
        self.anndata_ingest.generate_cluster_type_declaration(
            self.anndata_ingest.obtain_adata(), self.cluster_name
        )
        with open(self.cluster_filename) as header_file:
            header = header_file.readline().split("\t")
            self.assertEqual(
                ['TYPE', 'numeric', "numeric\n"],
                header,
                "did not find expected headers",
            )

    def test_generate_cluster_body(self):
        self.anndata_ingest.generate_cluster_body(
            self.anndata_ingest.obtain_adata(), self.cluster_name
        )
        compressed_file = self.cluster_filename + ".gz"
        with gzip.open(compressed_file, "rt", encoding="utf-8-sig") as cluster_body:
            line = cluster_body.readline().split("\t")
            expected_line = ['AAACATACAACCAC-1', '16.009954', "-21.073845\n"]
            self.assertEqual(
                expected_line,
                line,
                'did not get expected coordinates from cluster body',
            )

    def test_generate_metadata_file(self):
        self.anndata_ingest.generate_metadata_file(
            self.anndata_ingest.obtain_adata(), self.metadata_filename
        )
        compressed_file = self.metadata_filename + ".gz"
        with gzip.open(compressed_file, "rt", encoding="utf-8-sig") as metadata_body:
            name_line = metadata_body.readline().split("\t")
            expected_names = [
                'NAME',
                'n_genes',
                'percent_mito',
                'n_counts',
                'louvain',
                'donor_id',
                'biosample_id',
                'sex',
                'species',
                'species__ontology_label',
                'disease',
                'disease__ontology_label',
                'organ',
                'organ__ontology_label',
                'library_preparation_protocol',
                "library_preparation_protocol__ontology_label\n",
            ]
            self.assertEqual(
                expected_names, name_line, 'did not get expected headers from metadata body'
            )
            type_line = metadata_body.readline().split("\t")
            expected_types = [
                'TYPE',
                'NUMERIC',
                'NUMERIC',
                'NUMERIC',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                'GROUP',
                "GROUP\n",
            ]
            self.assertEqual(
                expected_types, type_line, 'did not get expected types from metadata body'
            )

    def test_gene_id_indexed_generate_processed_matrix(self):
        """Tests creating matrix when indexed by Ensembl ID, not gene name

        To reproduce test data:

        1.  Go to https://singlecell.broadinstitute.org/single_cell/study/SCP2557
        2.  In "Download" tab, click to download file "MRCA_AC.h5ad"
        3.  cp -p ~/Downloads/MRCA_AC.h5ad ~/scp-ingest-pipeline
        4.  source env/bin/activate
        5.  python
        6.  Run the following commands in your Python shell:

            import anndata
            adata = anndata.read_h5ad('MRCA_AC.h5ad')

            # Delete data not needed for this test
            del adata.obsm
            del adata.raw

            # Subset data to a small cohort of cells (males of age "P14")
            filtered_adata = adata[(adata.obs["sex"] == "male") & (adata.obs["age"] == "P14")]

            # Subset data to a small cohort of genes (those that are > 20 kbp in length)
            filtered_adata = filtered_adata[:, filtered_adata.var["feature_length"].astype('int64') > 20000]

            # Write out a file to disk for filtered AnnData
            filtered_adata.write('indexed_by_gene_id.h5ad')
        """
        indexed_by_geneid = AnnDataIngestor(
            "../tests/data/anndata/indexed_by_gene_id.h5ad", self.study_id, self.study_file_id
        )
        adata = indexed_by_geneid.obtain_adata()
        self.anndata_ingest.generate_processed_matrix(adata)

        now = time.time() # current time (ms since epoch)
        expected_features_fp = 'h5ad_frag.features.processed.tsv.gz'
        mtime = os.path.getmtime(expected_features_fp) # modified time (ms since epoch)
        self.assertTrue(abs(now - mtime) < 1000)

        with gzip.open(expected_features_fp, 'rt') as f:
            first_line = f.readline().strip('\n')

        expected_first_line = 'ENSMUSG00000037270\t4932438A13Rik'
        self.assertEqual(
            first_line, expected_first_line, 'Expected Ensembl ID and gene name'
        )


    def test_get_files_to_delocalize(self):
        files = AnnDataIngestor.clusterings_to_delocalize(self.valid_kwargs)
        compressed_file = self.cluster_filename + ".gz"
        expected_files = [compressed_file]
        self.assertEqual(expected_files, files)

    def test_delocalize_files(self):
        # just create header, no reason to run full extract
        self.anndata_ingest.generate_cluster_header(
            self.anndata_ingest.obtain_adata(), self.cluster_name
        )
        with patch('ingest_files.IngestFiles.delocalize_file'):
            AnnDataIngestor.delocalize_file(
                "gs://fake_bucket",
                self.study_id,
                AnnDataIngestor.clusterings_to_delocalize(self.valid_kwargs),
            )
            self.assertEqual(
                IngestFiles.delocalize_file.call_count,
                1,
                "expected 1 call to delocalize output files",
            )

