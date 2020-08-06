"""Integration tests for Genomes Pipeline; isolated tests for observable output

These tests verify that reference genome data can be extracted from upstream
sources (e.g. NCBI, Ensembl) transformed, as expected by code that loads
transformed data to Google Cloud Storage.

Test doubles are used for test speed and isolation.

PREREQUISITES
See https://github.com/broadinstitute/scp-ingest-pipeline#prerequisites

Also, `cd tests` before running the examples below.

EXAMPLES

# Run all tests
pytest

# Run all tests and see print() output
pytest -s

# Run only tests in test_genomes.py
pytest test_genomes.py

# Run all tests, show code coverage metrics
pytest --cov=../ingest/

"""

import unittest

# from unittest.mock import patch
import sys

sys.path.append('../ingest/genomes')
from genomes_pipeline import create_parser, parse_assemblies, parse_genome_annotations


class GenomesTestCase(unittest.TestCase):
    def setup_genomes(self, args):
        parsed_args = create_parser().parse_args(args)
        parse_assemblies(parsed_args)
        parse_genome_annotations(parsed_args)
        return

    def test_genomes_default(self):
        args = ['--input-dir', '../ingest/genomes/', '--local-output-dir', '.']
        self.setup_genomes(args)
