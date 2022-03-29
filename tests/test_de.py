""" test_de.py
    integration test to verify that de process generates expected output
"""

import unittest
import sys
import hashlib

sys.path.append("../ingest")
from cell_metadata import CellMetadata
from clusters import Clusters
from de import DifferentialExpression


class TestDifferentialExpression(unittest.TestCase):
    def test_de_process(self):
        """ Run DE on small test case
            confirm expected output
        """
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
            "method": "wilcoxon",
        }

        de = DifferentialExpression(
            cluster,
            cm,
            "../tests/data/differential_expression/de_integration.tsv",
            "dense",
            "cell_type__ontology_label",
            **de_kwargs,
        )
        de.execute_de()

        # md5 checksum calculated using reference file in tests/data/differential_expression/reference
        expected_checksum = "e3cc75eb3226ec8a2198205bc3e4581e"

        # running DifferentialExpression via pytest results in output files in the tests dir
        with open(
            "../tests/de_integration--cell_type__ontology_label--cholinergic_neuron--wilcoxon.tsv",
            "rb",
        ) as f:
            bytes = f.read()
            de_output_checksum = hashlib.md5(bytes).hexdigest()
        self.assertEqual(
            de_output_checksum,
            expected_checksum,
            "generated output file should match expected checksum",
        )

