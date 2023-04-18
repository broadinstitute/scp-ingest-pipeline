"""Integration tests for make_toy_data.py

EXAMPLES

See See https://github.com/broadinstitute/scp-ingest-pipeline#test.

"""

import unittest
import sys
import shutil
from random import randrange

sys.path.append('../ingest')
sys.path.append('../ingest/genomes')
from make_toy_data import create_parser, make_toy_data

output_dir = f'tmp/test_make_toy_{randrange(100_000)}/'


class MakeToyTestCase(unittest.TestCase):
    def execute_make_toy(self, args):
        """Matches main() in make_toy_data.py, but called as module"""
        parsed_args = create_parser().parse_args(args)
        make_toy_data(parsed_args)
        return

    @staticmethod
    def clean_up(output_dir):
        try:
            shutil.rmtree(output_dir)
        except OSError:
            print('no files to remove')

    def test_make_toy_data_default(self):
        args = ['--output-dir', output_dir]
        self.execute_make_toy(args)
        self.clean_up(output_dir)
