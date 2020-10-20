"""Integration tests for Genomes Pipeline

These tests verify that reference genome data can be extracted from upstream
sources (e.g. NCBI, Ensembl) transformed, as expected by code that loads
transformed data to Google Cloud Storage.

Test doubles are used for test speed and isolation.

EXAMPLES

See https://github.com/broadinstitute/scp-ingest-pipeline#test.
"""

import unittest

from unittest.mock import patch
import sys
import shutil
from random import randrange
from glob import glob

sys.path.append('../ingest/genomes')
from genomes_pipeline import create_parser, parse_assemblies, parse_genome_annotations


def mock_upload_ensembl_gtf_products(ensembl_metadata, scp_species, config):
    """Avoids call that requires credentials and long-running upload
    """
    return


def mock_record_annotation_metadata(output_dir, ensembl_metadata, scp_species):
    """Avoids call that (upstream) requires credentials and long-running upload

    To consider:
    Refactor `record_annotation_metadata` to be *before* upload_ensembl_gtf_products,
    then remove this mock.
    """
    return


output_dir = f'tmp/test_genomes_{randrange(100_000)}/'


class GenomesTestCase(unittest.TestCase):
    @patch(
        'genome_annotation_metadata.upload_ensembl_gtf_products',
        side_effect=mock_upload_ensembl_gtf_products,
    )
    @patch(
        'genome_annotation_metadata.record_annotation_metadata',
        side_effect=mock_record_annotation_metadata,
    )
    def execute_genomes_pipeline(
        self, args, mock_record_annotation_metadata, mock_upload_ensembl_gtf_products
    ):
        """Matches main() in genomes_pipeline.py, but called as module
        """
        parsed_args = create_parser().parse_args(args)
        parse_assemblies(parsed_args)
        parse_genome_annotations(parsed_args)
        return

    @staticmethod
    def clean_up(output_dir):
        try:
            shutil.rmtree(output_dir)
        except OSError:
            print('no files to remove')

    def test_genomes_default(self):
        """Genomes Pipeline should extract and transform reference data
        """
        args = ['--input-dir', 'mock_data/genomes/', '--local-output-dir', output_dir]
        self.execute_genomes_pipeline(args)

        # There should be two gzipped GTF files
        gzipped_gtfs = glob(f'{output_dir}*.gtf.gz')
        self.assertEqual(len(gzipped_gtfs), 2)

        self.clean_up(output_dir)
