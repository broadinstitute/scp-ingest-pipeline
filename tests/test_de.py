""" test_de.py
    integration test to verify that de process generates expected output
"""

import unittest
import sys
import hashlib
import os
import pandas as pd

sys.path.append("../ingest")
from cell_metadata import CellMetadata
from clusters import Clusters
from de import DifferentialExpression


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


def find_expected_files(labels, cluster_name, annotation, method):
    """ Check that files were created for all expected annotation labels
    """
    found = 0
    for label in labels:
        sanitized_label = label.replace(" ", "_")
        expected_file = (
            f"{cluster_name}--{annotation}--{str(sanitized_label)}--{method}.tsv"
        )
        assert os.path.exists(expected_file)
        found += 1
    return found


class TestDifferentialExpression(unittest.TestCase):
    def test_process_missing_metadata(self):
        cm = CellMetadata(
            "../tests/data/differential_expression/de_integration_unordered_metadata.tsv",
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
            "../tests/data/differential_expression/de_integration_unordered_metadata.tsv",
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

    def test_de_process_dense(self):
        """ Run DE on small test case with dense matrix inputs
            confirm expected output
        """
        test_annotation = "cell_type__ontology_label"
        test_method = "wilcoxon"
        cm = CellMetadata(
            "../tests/data/differential_expression/de_integration_unordered_metadata.tsv",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            study_accession="SCPde",
            tracer=None,
        )

        cluster = Clusters(
            "../tests/data/differential_expression/de_integration_cluster.tsv",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            "de_integration",
        )

        de_kwargs = {
            "study_accession": cm.study_accession,
            "name": cluster.name,
            "method": test_method,
        }

        de = DifferentialExpression(
            cluster,
            cm,
            "../tests/data/differential_expression/de_integration.tsv",
            "dense",
            test_annotation,
            **de_kwargs,
        )
        de.execute_de()
        de_cells = DifferentialExpression.get_cluster_cells(cluster.file['NAME'].values)
        labels = get_annotation_labels(cm, test_annotation, de_cells)
        found_label_count = find_expected_files(
            labels, cluster.name, test_annotation, test_method
        )

        self.assertEqual(
            found_label_count,
            5,
            f"expected five annotation labels for {test_annotation}",
        )

        expected_file_path = (
            "../tests/de_integration--cell_type__ontology_label"
            "--cholinergic_neuron--wilcoxon.tsv"
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

    def test_de_process_sparse(self):
        """ Run DE on small test case with sparse matrix inputs
                confirm expected output
            """
        test_annotation = "cell_type__ontology_label"
        test_method = "wilcoxon"
        cm = CellMetadata(
            "../tests/data/differential_expression/sparse/sparsemini_metadata.txt",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            study_accession="SCPsparsemini",
            tracer=None,
        )

        cluster = Clusters(
            "../tests/data/differential_expression/sparse/sparsemini_cluster.txt",
            "addedfeed000000000000000",
            "dec0dedfeed0000000000000",
            "de_sparse_integration",
        )

        de_kwargs = {
            "study_accession": cm.study_accession,
            "name": cluster.name,
            "method": test_method,
            "gene_file": "../tests/data/differential_expression/sparse/sparsemini_features.tsv",
            "barcode_file": "../tests/data/differential_expression/sparse/sparsemini_barcodes.tsv",
        }

        de = DifferentialExpression(
            cluster,
            cm,
            "../tests/data/differential_expression/sparse/sparsemini_matrix.mtx",
            "mtx",
            test_annotation,
            **de_kwargs,
        )
        de.execute_de()
        de_cells = DifferentialExpression.get_cluster_cells(cluster.file['NAME'].values)
        labels = get_annotation_labels(cm, test_annotation, de_cells)
        # In find_expected_files, checks all files with expected names were created
        # yields the number of files expected for an external check for file count
        found_label_count = find_expected_files(
            labels, cluster.name, test_annotation, test_method
        )

        self.assertEqual(
            found_label_count,
            7,
            f"expected seven annotation labels for {test_annotation}",
        )

        expected_file_path = (
            "../tests/de_sparse_integration--cell_type__ontology_label"
            "--endothelial_cell--wilcoxon.tsv"
        )

        content = pd.read_csv(expected_file_path, sep="\t", index_col=0)
        # confirm expected gene in DE file at expected position
        self.assertEqual(
            content.iloc[1, 0],
            "Mrpl15",
            "Did not find expected gene, Mrpl15, at second row in DE file",
        )
        # confirm calculated value has expected significant digits
        self.assertEqual(
            content.iloc[0, 2],
            11.63,
            "Did not find expected logfoldchange value for Sox17 in DE file",
        )

        # md5 checksum calculated using reference file in tests/data/differential_expression/sparse/reference
        expected_checksum = "07b6c6565430a17f4f048e7b4f53ddac"

        # running DifferentialExpression via pytest results in output files in the tests dir
        with open(expected_file_path, "rb") as f:
            bytes = f.read()
            de_output_checksum = hashlib.md5(bytes).hexdigest()
        self.assertEqual(
            de_output_checksum,
            expected_checksum,
            "generated output file should match expected checksum",
        )

