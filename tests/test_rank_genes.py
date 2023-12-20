"""Test author_de.py

# Run all tests in a manner that shows report_issues output
pytest test_author_de.py -s
"""

import unittest
import sys
import glob
import os
import pytest
from unittest.mock import patch
from unittest.mock import Mock
from requests import Response

sys.path.append("../ingest")
from rank_genes import fetch_publication_text, fetch_gene_cache, fetch_context

mock_response = Mock(spec=Response)

class TestRankGenes(unittest.TestCase):

    # def test_publication_text(self):
    #     publication_text = fetch_publication_text("https://doi.org/10.1126/sciadv.adf6251")
    #     self.assertEqual(len(publication_text), 66433, "Publication text does not have expected length")

    # def test_fetch_gene_cache(self):
    #     gene_cache = fetch_gene_cache('homo-sapiens')
    #     interest_rank_by_gene = gene_cache[0]
    #     error_msg = "Human should have > 20K interest-ranked genes"
    #     self.assertGreater(len(interest_rank_by_gene.keys()), 20000, error_msg)

    def test_fetch_context(self):
        expected_json = {
            "bucketId": "fc-1234",
            "taxonNames": ["Homo sapiens"]
        }
        mock_response.json = lambda: {"ok": "true"}
        with patch("requests.get", return_value=mock_response):
            bucket_id, organism, de_dict = fetch_context("SCP123")
            self.assertEqual(bucket_id, 'adsf')



    def teardown_method(self, test_method):
        files = glob.glob('cluster_umap_txt--General_Celltype*.tsv')
        for file in files:
            os.remove(file)



