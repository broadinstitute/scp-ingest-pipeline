"""Integration tests for make_toy_data.py

PREREQUISITES
See https://github.com/broadinstitute/scp-ingest-pipeline#prerequisites

Also, `cd tests` before running the examples below.

EXAMPLES

# Run all tests
pytest

# Run all tests and see print() output
pytest -s

# Run only tests in test_genomes.py
pytest test_make_toy.py

# Run all tests, show code coverage metrics
pytest --cov=../ingest/

"""

import unittest
import sys

sys.path.append('../ingest')
sys.path.append('../ingest/genomes')
from make_toy_data import create_parser, make_toy_data


class MakeToyTestCase(unittest.TestCase):
    def execute_make_toy(self, args):
        parsed_args = create_parser().parse_args(args)
        make_toy_data(parsed_args)
        return

    def test_make_toy_data_default(self):
        args = []
        self.execute_make_toy(args)
