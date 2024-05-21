"""Test author_de.py

# Run all tests in a manner that shows report_issues output
pytest test_author_de.py -s
"""

import unittest
import sys
import glob
import os
import pytest
import json
from unittest.mock import patch
from unittest.mock import Mock
from requests import Response

sys.path.append("../ingest")
from rank_genes import fetch_publication_text, fetch_gene_cache, fetch_context

mock_response = Mock(spec=Response)

class TestRankGenes(unittest.TestCase):

    def test_publication_text(self):
        [publication_text, platform] = fetch_publication_text("https://doi.org/10.1126/sciadv.adf6251")
        self.assertEqual(len(publication_text), 66433, "Publication text does not have expected length")

    def test_fetch_gene_cache(self):
        gene_cache = fetch_gene_cache('homo-sapiens')
        interest_rank_by_gene = gene_cache[0]
        error_msg = "Human should have > 20K interest-ranked genes"
        self.assertGreater(len(interest_rank_by_gene.keys()), 20000, error_msg)

    def test_fetch_context(self):
        with open("../tests/data/rank_genes/mock_explore_response.json") as f:
            expected_json = json.loads(f.read())
        mock_response.status_code = 200
        mock_response.json = lambda: expected_json
        with patch(
            "requests.get", return_value=mock_response
        ), patch.dict(
            "os.environ", {"DATABASE_NAME": "foo_production"}
        ):
            bucket_id, organism, de_dict = fetch_context("SCP123")
            self.assertEqual(bucket_id, "fc-1a23b456-cde4-5fa6-bc7d-89e0f1a23bc4")
            self.assertEqual(organism, 'homo-sapiens')
            self.assertEqual(de_dict["de_groups_and_files"][3][0], "Fibroblasts")

    def teardown_method(self, test_method):
        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        for file in files:
            os.remove(file)



