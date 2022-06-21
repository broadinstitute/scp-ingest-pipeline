""" test_de.py
    integration test to verify that de process generates expected output
"""

import unittest
import sys
import hashlib
import os
import glob
import pandas as pd
from unittest.mock import patch
import scanpy as sc
import re

sys.path.append("../ingest")
from cell_metadata import CellMetadata
from clusters import Clusters
from de import DifferentialExpression
from ingest_files import IngestFiles


def get_annotation_labels(metadata, annotation, de_cells):
    """ extract annotation_labels for specified annotation
        from metadata file filtered by cluster file cells
    """
    dtypes = DifferentialExpression.determine_dtypes(
        metadata.headers, metadata.annot_types
    )
    de_annots = DifferentialExpression.subset_annots(metadata, de_cells)
    unique_labels = de_annots[annotation].unique()
    return unique_labels.tolist()


def find_expected_files(labels, cluster_name, annotation, scope, method):
    """ Check that files were created for all expected annotation labels
    """
    found = []
    sanitized_cluster_name = re.sub(r'\W', '_', cluster_name)
    sanitized_annotation = re.sub(r'\W', '_', annotation)
    for label in labels:
        sanitized_label = re.sub(r'\W', '_', label)
        expected_file = f"{sanitized_cluster_name}--{sanitized_annotation}--{sanitized_label}--{scope}--{method}.tsv"
        assert os.path.exists(expected_file)
        found.append(expected_file)
    return found


def run_de(**test_config):
    test_annotation = test_config["test_annotation"]
    test_scope = test_config["test_scope"]
    test_method = test_config["test_method"]
    annot_path = test_config["annot_path"]
    study_accession = test_config["study_accession"]
    cluster_path = test_config["cluster_path"]
    cluster_name = test_config["cluster_name"]
    matrix_path = test_config["matrix_file"]
    matrix_type = test_config["matrix_type"]

    cm = CellMetadata(
        annot_path,
        "addedfeed000000000000000",
        "dec0dedfeed0000000000000",
        study_accession=study_accession,
        tracer=None,
    )

    cluster = Clusters(
        cluster_path,
        "addedfeed000000000000000",
        "dec0dedfeed0000000000000",
        cluster_name,
    )

    de_kwargs = {
        "study_accession": cm.study_accession,
        "name": cluster.name,
        "annotation_scope": test_scope,
        "method": test_method,
    }

    if "gene_file" in test_config:
        de_kwargs["gene_file"] = test_config["gene_file"]

    if "barcode_file" in test_config:
        de_kwargs["barcode_file"] = test_config["barcode_file"]

    de = DifferentialExpression(
        cluster, cm, matrix_path, matrix_type, test_annotation, **de_kwargs
    )
    de.execute_de()
    de_cells = DifferentialExpression.get_cluster_cells(cluster.file['NAME'].values)
    labels = get_annotation_labels(cm, test_annotation, de_cells)
    # In find_expected_files, checks all files with expected names were created
    # yields the number of files expected for an external check for file count
    found_labels = find_expected_files(
        labels, cluster.name, test_annotation, test_scope, test_method
    )
    return found_labels


class TestDifferentialExpression(unittest.TestCase):
    def test_process_missing_metadata(self):
        cm = CellMetadata(
            "../tests/data/differential_expression/de_dense_metadata.tsv",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            study_accession="SCPde",
            tracer=None,
        )
        dtypes = DifferentialExpression.determine_dtypes(cm.headers, cm.annot_types)
        annots = DifferentialExpression.process_annots(
            cm.file_path, cm.ALLOWED_FILE_TYPES, cm.headers, dtypes
        )
        annot_labels = annots["cell_type__ontology_label"].unique().tolist()
        self.assertIn(
            "__Unspecified__",
            annot_labels,
            "Expected NaN conversion to \"__Unspecified__\"",
        )

    def test_assess_annotation(self):
        test_annotation = "seurat_clusters"
        cm = CellMetadata(
            "../tests/data/differential_expression/de_dense_metadata.tsv",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            study_accession="SCPde",
            tracer=None,
        )
        # alter metadata so seurat_clusters is TYPE numeric
        cm.annot_types[21] = 'numeric'
        self.assertRaises(
            TypeError, DifferentialExpression.assess_annotation, test_annotation, cm
        )
        # alter metadata so seurat_clusters has bad TYPE designation
        cm.annot_types[21] = 'foo'
        self.assertRaises(
            ValueError, DifferentialExpression.assess_annotation, test_annotation, cm
        )
        # remove seurat_clusters from metadata
        cm.headers.pop(21)
        self.assertRaises(
            KeyError, DifferentialExpression.assess_annotation, test_annotation, cm
        )

    def test_detect_duplicate_gene_names(self):
        """Genes file can have one or two columns of gene information
            If two columns present, use the second column containing gene names
            unless there are duplicate gene names in the second column
            If duplicates, check that 1st plus 2nd column provides uniqueness
            If unique when joined, join columns with pipe (|) for use as DE input
        """
        no_dup_genes_path = (
            "../tests/data/differential_expression/sparse/sparsemini_features.tsv"
        )
        no_dup_genes = DifferentialExpression.get_genes(no_dup_genes_path)
        self.assertNotIn(
            "|", no_dup_genes[0], f"no delimiter expected in {no_dup_genes[0]}"
        )

        dup_genes_path = (
            "../tests/data/differential_expression/sparse/sparsemini_dup_gene_name.tsv"
        )
        dup_genes = DifferentialExpression.get_genes(dup_genes_path)
        self.assertIn("|", dup_genes[0], f"no delimiter expected in {dup_genes[0]}")

    def test_delimiter_in_gene_name(self):
        delimited_data = {"names": ["Tns1", "Gfra1"], "scores": ["10.5", "10.34"]}
        delimited_df = pd.DataFrame(delimited_data)
        self.assertFalse(
            DifferentialExpression.delimiter_in_gene_name(delimited_df),
            "no pipe delimiter should be detected in the input",
        )

        undelimited_data = {
            "names": ["ENSMUST00000027035|Sox17", "ENSMUST00000195555|Sox17"],
            "scores": ["41.459137", "-5.058518"],
        }
        undelimited_df = pd.DataFrame(undelimited_data)
        self.assertTrue(
            DifferentialExpression.delimiter_in_gene_name(undelimited_df),
            "expected pipe delimiter undetected",
        )

    def test_filename_sanitation(self):
        """ Bugfix (SCP-4459) so sanitization does not collapse adjacent non-alphanumeric characters to
            single underscores, see also SCP-4455 for manual fix
        """
        arguments = {
            "cluster_name": "UMAP, pre-QC all cells (complexity greater than or equal to 1000)",
            "annotation_name": "cell..type",
        }
        files_to_match = DifferentialExpression.string_for_output_match(arguments)
        print(files_to_match)
        self.assertEqual(
            files_to_match,
            "UMAP__pre_QC_all_cells__complexity_greater_than_or_equal_to_1000_--cell__type*.tsv",
            "unexpected result from sanitation function",
        )

    def test_de_remove_single_sample(self):
        """ Test single sample removal
        """
        test_annotation = "seurat_clusters"
        cm = CellMetadata(
            "../tests/data/differential_expression/de_singlesample_metadata.tsv",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            study_accession="SCPde",
            tracer=None,
        )
        adata = sc.read("../tests/data/differential_expression/de_dense_matrix.tsv")
        adata = adata.transpose()
        dtypes = DifferentialExpression.determine_dtypes(cm.headers, cm.annot_types)
        annots = DifferentialExpression.process_annots(
            cm.file_path, cm.ALLOWED_FILE_TYPES, cm.headers, dtypes
        )
        adata.obs = DifferentialExpression.order_annots(annots, adata.obs_names)
        orig_obs = adata.n_obs
        adata2 = DifferentialExpression.remove_single_sample_data(
            adata, test_annotation
        )
        new_obs = adata2.n_obs
        number_removed = orig_obs - new_obs
        self.assertEqual(
            number_removed,
            1,
            f"expected removal of single sample, found {number_removed}",
        )

        removed = list(set(adata.obs_names) - set(adata2.obs_names))
        self.assertEqual(
            removed,
            ["Mm_AMB_N107"],
            f"expected removal of cell [\'Mm_AMB_N107\'], not \"{removed}\"",
        )

    def test_de_process_dense(self):
        """ Run DE on small test case with dense matrix inputs
            confirm expected output
        """
        test_annotation = "cell_type__ontology_label"
        test_config = {
            "test_annotation": test_annotation,
            "test_scope": "study",
            "test_method": "wilcoxon",
            "annot_path": "../tests/data/differential_expression/de_dense_metadata.tsv",
            "study_accession": "SCPde",
            "cluster_path": "../tests/data/differential_expression/de_dense_cluster.tsv",
            "cluster_name": "de_integration",
            "matrix_file": "../tests/data/differential_expression/de_dense_matrix.tsv",
            "matrix_type": "dense",
        }

        found_labels = run_de(**test_config)
        found_label_count = len(found_labels)

        self.assertEqual(
            found_label_count,
            5,
            f"expected five annotation labels for {test_annotation}",
        )

        expected_file = (
            "de_integration--cell_type__ontology_label"
            "--cholinergic_neuron--study--wilcoxon.tsv"
        )

        expected_file_path = f"../tests/{expected_file}"

        # confirm expected results filename was generated in found result files
        self.assertIn(
            expected_file, found_labels, "Expected filename not in found files list"
        )

        content = pd.read_csv(expected_file_path, sep="\t", index_col=0)
        # confirm expected gene in DE file at expected position
        self.assertEqual(
            content.iloc[3, 0],
            "Vipr2",
            "Did not find expected gene, Vipr2, at fourth row in DE file",
        )
        # confirm calculated value has expected significant digits
        self.assertEqual(
            content.iloc[183, 2],
            -4.834,
            "Did not find expected logfoldchange value for Nsg2 in DE file",
        )

        # md5 checksum calculated using reference file in tests/data/differential_expression/reference
        expected_checksum = "e3cc75eb3226ec8a2198205bc3e4581e"

        # running DifferentialExpression via pytest results in output files in the tests dir
        with open(expected_file_path, "rb") as f:
            bytes = f.read()
            de_output_checksum = hashlib.md5(bytes).hexdigest()
        self.assertEqual(
            de_output_checksum,
            expected_checksum,
            "generated output file should match expected checksum",
        )

        # clean up DE outputs
        output_wildcard_match = f"../tests/de_integration--{test_annotation}*.tsv"
        files = glob.glob(output_wildcard_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_de_process_sparse(self):
        """ Run DE on small test case with sparse matrix inputs
                confirm expected output
            """
        test_annotation = "cell_type__ontology_label"
        test_config = {
            "test_annotation": test_annotation,
            "test_scope": "study",
            "test_method": "wilcoxon",
            "annot_path": "../tests/data/differential_expression/sparse/sparsemini_metadata.txt",
            "study_accession": "SCPsparsemini",
            "cluster_path": "../tests/data/differential_expression/sparse/sparsemini_cluster.txt",
            "cluster_name": "de_sparse_dup_gene",
            "matrix_file": "../tests/data/differential_expression/sparse/sparsemini_matrix.mtx",
            "matrix_type": "mtx",
            "gene_file": "../tests/data/differential_expression/sparse/sparsemini_dup_gene_name.tsv",
            "barcode_file": "../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv",
        }

        found_labels = run_de(**test_config)
        found_label_count = len(found_labels)

        self.assertEqual(
            found_label_count,
            7,
            f"expected seven annotation labels for {test_annotation}",
        )

        expected_file = (
            "de_sparse_dup_gene--cell_type__ontology_label"
            "--endothelial_cell--study--wilcoxon.tsv"
        )

        expected_file_path = f"../tests/{expected_file}"

        # confirm expected results filename was generated in found result files
        self.assertIn(
            expected_file, found_labels, "Expected filename not in found files list"
        )

        content = pd.read_csv(expected_file_path, sep="\t", index_col=0)
        # confirm expected gene in DE file at expected position
        self.assertEqual(
            content.iloc[1, 0],
            "Sox17",
            "Did not find expected gene, Sox17, at second row in DE file.",
        )
        # confirm calculated value has expected significant digits
        self.assertEqual(
            content.iloc[0, 2],
            11.63,
            "Did not find expected logfoldchange value for Sox17 in DE file.",
        )
        # confirm duplicate gene input generates expected gene_id info in output
        self.assertIn('gene_id', content.columns, "Expected gene_id output not found.")

        # md5 checksum calculated using reference file in tests/data/differential_expression/sparse/reference
        # file updated 2022-05-25 to include output for duplicate gene handling
        expected_checksum = "ca0c7dcc4048614f22d6bc7dec18a2c0"

        # running DifferentialExpression via pytest results in output files in the tests dir
        with open(expected_file_path, "rb") as f:
            bytes = f.read()
            de_output_checksum = hashlib.md5(bytes).hexdigest()
        self.assertEqual(
            de_output_checksum,
            expected_checksum,
            "Generated output file should match expected checksum.",
        )

        expected_output_match = "de_sparse_dup_gene--cell_type__ontology_label*.tsv"

        with patch('ingest_files.IngestFiles.delocalize_file'):
            DifferentialExpression.delocalize_de_files(
                'gs://fake_bucket', None, expected_output_match
            )

            self.assertEqual(
                IngestFiles.delocalize_file.call_count,
                7,
                "expected 7 calls to delocalize output files",
            )

        # clean up DE outputs
        files = glob.glob(expected_output_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_de_process_na(self):
        """ Run DE on small test case with na-type values in matrix
            confirm expected output filenames
        """
        test_annotation = "cell_type__ontology_label"
        test_config = {
            "test_annotation": test_annotation,
            "test_scope": "study",
            "test_method": "wilcoxon",
            "annot_path": "../tests/data/differential_expression/de_dense_metadata_na.txt",
            "study_accession": "SCPna",
            "cluster_path": "../tests/data/differential_expression/de_dense_cluster.tsv",
            "cluster_name": "de_na",
            "matrix_file": "../tests/data/differential_expression/de_dense_matrix.tsv",
            "matrix_type": "dense",
        }

        found_labels = run_de(**test_config)
        found_label_count = len(found_labels)

        self.assertEqual(
            found_label_count,
            9,
            f"expected nine annotation labels for {test_annotation}",
        )

        expected_file = "de_na--cell_type__ontology_label--N_A--study--wilcoxon.tsv"

        # confirm expected results filename was generated in found result files
        self.assertIn(
            expected_file, found_labels, "Expected filename not in found files list"
        )

        expected_output_match = "de_na--cell_type__ontology_label*.tsv"

        # clean up DE outputs
        files = glob.glob(expected_output_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")

    def test_de_process_sanitize(self):
        """ Run DE on small test case with na-type values in matrix
            confirm expected output filenames
        """
        test_annotation = "misc++cellaneous"
        test_config = {
            "test_annotation": test_annotation,
            "test_scope": "study",
            "test_method": "wilcoxon",
            "annot_path": "../tests/data/differential_expression/de_dense_metadata_sanitize.txt",
            "study_accession": "SCPsanitize",
            "cluster_path": "../tests/data/differential_expression/de_dense_cluster.tsv",
            "cluster_name": "UMAP, pre-QC",
            "matrix_file": "../tests/data/differential_expression/de_dense_matrix.tsv",
            "matrix_type": "dense",
        }

        found_labels = run_de(**test_config)
        found_label_count = len(found_labels)

        self.assertEqual(
            found_label_count,
            5,
            f"expected five annotation labels for {test_annotation}",
        )

        expected_file = (
            "UMAP__pre_QC--misc__cellaneous--cholinergic__neuron_--study--wilcoxon.tsv"
        )

        # confirm expected results filename was generated in found result files
        self.assertIn(
            expected_file, found_labels, "Expected filename not in found files list"
        )

        expected_output_match = "UMAP__pre_QC--misc__cellaneous*.tsv"

        # clean up DE outputs
        files = glob.glob(expected_output_match)

        for file in files:
            try:
                os.remove(file)
            except:
                print(f"Error while deleting file : {file}")
