"""Test validation/minify_ontology.py

# Run all tests in a manner that shows report_issues output
pytest test_minify_ontology.py -s
"""

import unittest
import sys
import glob
import os
import gzip
from unittest.mock import patch
from unittest.mock import Mock

sys.path.append("../ingest/validation")
from minify_ontology import OntologyMinifier

# mock_response = Mock(spec=Response)

class TestOntologyMinifier(unittest.TestCase):

    def test_mondo_and_pato_minification(self):
        OntologyMinifier(['disease'], False)
        files = glob.glob('*.tsv.gz')
        self.assertEqual(len(files), 2, 'Did not find 2 TSV.GZ files')
        with gzip.open('mondo.min.tsv.gz', 'rt') as f:
            first_line = f.readline().strip().split('\t')
        expected_first_line = [
            'MONDO_0000001',
            'disease',
            'condition||disease||disease or disorder||disease or disorder, non-neoplastic||diseases||diseases and disorders||disorder||disorders||medical condition||other disease'
        ]
        error_message = 'Did not get expected first line in mondo.min.tsv.gz'
        self.assertEqual(first_line, expected_first_line, error_message)

    def teardown_method(self, test_method):
        output_files = glob.glob('*.min.tsv.gz')
        for file in output_files:
            os.remove(file)



