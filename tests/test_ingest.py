"""Integration tests for Ingest Pipeline; isolated tests for observable output

These tests verify that various matrix file types can be extracted and
transformed, as expected by code that loads transformed data into MongoDB.

Test doubles are used for test speed and isolation.

PREREQUISITE
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
Also, before 'EXAMPLES' commands, run `cd tests`.

EXAMPLES

# Run all tests
pytest

# Run all tests and see print() output
pytest -s

# Run only tests in test_ingest.py
pytest test_ingest.py

# Run only the test "test_good_metadata_file" in this test module
pytest test_ingest.py -k 'test_good_metadata_file' -s

# Run all tests, using multiple CPUs
# Details: https://pypi.org/project/pytest-xdist/
pytest -n auto -s

# Run all tests, show code coverage metrics
pytest --cov=../ingest/

"""
import ast
import sys
import unittest
from unittest.mock import patch
from bson.objectid import ObjectId
from mock_data.dense_matrix_19_genes_100k_cells_txt.gene_models_0 import gene_models
from mock_data.matrix_mtx.gene_model_0 import expected_model
from gcp_mocks import mock_storage_client, mock_storage_blob

sys.path.append('../ingest')
from ingest_pipeline import (
    create_parser,
    validate_arguments,
    IngestPipeline,
    exit_pipeline,
    run_ingest,
)


def mock_load(self, *args, **kwargs):
    """Enables overwriting normal function with this placeholder.
    Returning the arguments enables tests to verify that the code invokes
    this method with expected argument values.

    TODO:
    Integrate MongoDB emulator for faster, higher-coverage tests (SCP-2000)

    This will enable us to also verify (and thus cover) loading-code *outputs*,
    unlike here where we merely give a way to verify loading-code *inputs*.
    Doing so via integration tests will isolate us from implementation changes.
    """
    self.load_args = args
    self.load_kwargs = kwargs


# Mock method that writes to database
IngestPipeline.load_expression_file = mock_load
IngestPipeline.load = mock_load


def get_gene_model(mock_dir):
    """Return actual and expected gene model, using actual and mock data
    """

    with open(f'mock_data/{mock_dir}/gene_model_0.txt') as f:
        # Create a dictionary from the string-literal mock
        expected_model = ast.literal_eval(f.read())
        # convert strings to BSON ObjectIds
        study_id = ObjectId(expected_model['study_id'])
        study_file_id = ObjectId(expected_model['study_file_id'])
        expected_model['study_id'] = study_id
        expected_model['study_file_id'] = study_file_id

    return expected_model


class IngestTestCase(unittest.TestCase):
    @patch('google.cloud.storage.Blob', side_effect=mock_storage_blob)
    @patch('google.cloud.storage.Client', side_effect=mock_storage_client)
    def setup_ingest(self, args, mock_storage_client, mock_storage_blob):

        self.maxDiff = None

        parsed_args = create_parser().parse_args(args)
        validate_arguments(parsed_args)
        arguments = vars(parsed_args)

        ingest = IngestPipeline(**arguments)

        status, status_cell_metadata = run_ingest(ingest, arguments, parsed_args)

        print(f'status is {status}')

        return ingest, arguments, status, status_cell_metadata

    # def test_ingest_dense_matrix(self):
    #     """Ingest Pipeline should extract, transform, and load dense matrices
    #     """
    #
    #     args = [
    #         '--study-id',
    #         '5d276a50421aa9117c982845',
    #         '--study-file-id',
    #         '5dd5ae25421aa910a723a337',
    #         'ingest_expression',
    #         '--taxon-name',
    #         'Homo sapiens',
    #         '--taxon-common-name',
    #         'human',
    #         '--ncbi-taxid',
    #         '9606',
    #         '--genome-assembly-accession',
    #         'GCA_000001405.15',
    #         '--genome-annotation',
    #         'Ensembl 94',
    #         '--matrix-file',
    #         'gs://fake-bucket/tests/data/dense_matrix_19_genes_1000_cells.txt',
    #         '--matrix-file-type',
    #         'dense',
    #     ]
    #     ingest = self.setup_ingest(args)[0]
    #     models = ingest.load_args[0]
    #     print(models)
    #     for model in models:
    #         # Ensure that 'ObjectID' in model is removed
    #         del model['_id']
    #         self.assertEqual(model, gene_models[model['name']])
    #     # print(models)
    #
    #     # Verify gene model looks as expected
    #     # mock_dir = 'dense_matrix_19_genes_100k_cells_txt'
    #     # expected_model = get_gene_model(mock_dir)
    #
    #     # self.assertEqual(models, gene_models)
    #
    # def test_ingest_local_dense_matrix(self):
    #     """Ingest Pipeline should extract and transform local dense matrices
    #     """
    #
    #     args = [
    #         '--study-id',
    #         '5d276a50421aa9117c982845',
    #         '--study-file-id',
    #         '5dd5ae25421aa910a723a337',
    #         'ingest_expression',
    #         '--taxon-name',
    #         'Homo sapiens',
    #         '--taxon-common-name',
    #         'human',
    #         '--ncbi-taxid',
    #         '9606',
    #         '--genome-assembly-accession',
    #         'GCA_000001405.15',
    #         '--genome-annotation',
    #         'Ensembl 94',
    #         '--matrix-file',
    #         '../tests/data/dense_matrix_19_genes_1000_cells.txt',
    #         '--matrix-file-type',
    #         'dense',
    #     ]
    #     ingest = self.setup_ingest(args)[0]
    #
    #     models = ingest.load_args[0]
    #     for model in models:
    #         # Ensure that 'ObjectID' in model is removed
    #         del model['_id']
    #         self.assertEqual(model, gene_models[model['name']])
    #     # print(models)
    #
    # def test_ingest_local_compressed_dense_matrix(self):
    #     """Ingest Pipeline should extract and transform local dense matrices
    #         from compressed file in the same manner as uncompressed file
    #     """
    #
    #     args = [
    #         '--study-id',
    #         '5d276a50421aa9117c982845',
    #         '--study-file-id',
    #         '5dd5ae25421aa910a723a337',
    #         'ingest_expression',
    #         '--taxon-name',
    #         'Homo sapiens',
    #         '--taxon-common-name',
    #         'human',
    #         '--ncbi-taxid',
    #         '9606',
    #         '--genome-assembly-accession',
    #         'GCA_000001405.15',
    #         '--genome-annotation',
    #         'Ensembl 94',
    #         '--matrix-file',
    #         '../tests/data/dense_matrix_19_genes_100k_cells.txt.gz',
    #         '--matrix-file-type',
    #         'dense',
    #     ]
    #     ingest = self.setup_ingest(args)[0]
    #
    #     models = ingest.load_args[0]
    #     for model in models:
    #         # Ensure that 'ObjectID' in model is removed
    #         del model['_id']
    #         self.assertEqual(model, gene_models[model['name']])

    def test_ingest_mtx_matrix(self):
        """Ingest Pipeline should extract and transform MTX matrix bundles
        """

        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/AB_toy_data_toy.matrix.mtx',
            '--matrix-file-type',
            'mtx',
            '--gene-file',
            '../tests/data/AB_toy_data_toy.genes.tsv',
            '--barcode-file',
            '../tests/data/AB_toy_data_toy.barcodes.tsv',
        ]
        ingest = self.setup_ingest(args)[0]
        models = ingest.load_args[0]
        for model in models:
            # Ensure that 'ObjectID' in model is removed
            del model['_id']
        # print(model)
        self.assertEqual(models, expected_model)

    def test_remote_mtx_bundles(self):
        """Ingest Pipeline should handle MTX matrix files fetched from bucket
        """

        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            'gs://fake-bucket/tests/data/AB_toy_data_toy.matrix.mtx',
            '--matrix-file-type',
            'mtx',
            '--gene-file',
            'gs://fake-bucket/tests/data/AB_toy_data_toy.genes.tsv',
            '--barcode-file',
            'gs://fake-bucket/tests/data/AB_toy_data_toy.barcodes.tsv',
        ]
        ingest, arguments, status, status_cell_metadata =self.setup_ingest(args)

        models = ingest.load_args[0]
        print(models)
        for model in models:
            # Ensure that 'ObjectID' in model is removed
            del model['_id']
            print(model)
        # print(model)
        self.assertEqual(models, expected_model)

    def test_mtx_bundle_argument_validation(self):
        """Omitting --gene-file and --barcode-file in MTX ingest should error
        """

        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_expression',
            '--taxon-name',
            'Homo sapiens',
            '--taxon-common-name',
            'human',
            '--ncbi-taxid',
            '9606',
            '--genome-assembly-accession',
            'GCA_000001405.15',
            '--genome-annotation',
            'Ensembl 94',
            '--matrix-file',
            '../tests/data/matrix.mtx',
            '--matrix-file-type',
            'mtx',
        ]

        self.assertRaises(ValueError, self.setup_ingest, args)

        # TODO: This test does not run.  De-indent and fix.
        def test_bad_format_dense(self):
            args = [
                '--study-id',
                '5d276a50421aa9117c982845',
                '--study-file-id',
                '5dd5ae25421aa910a723a337',
                'ingest_expression',
                '--matrix-file',
                '../tests/data/expression_matrix_bad_missing_keyword.txt',
                '--matrix-file-type',
                'dense',
            ]
            with self.assertRaises(SystemExit) as cm:
                self.setup_ingest(args)
            not self.assertEqual(cm.exception.code, 0)

    def test_good_metadata_file(self):
        """Ingest Pipeline should succeed for properly formatted metadata file
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cell_metadata',
            '--cell-metadata-file',
            '../tests/data/metadata_example.txt',
            '--study-accession',
            'SCP123',
            '--ingest-cell-metadata',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        # TODO:
        # After integrating MongoDB emulator (SCP-2000), refactor this test to
        # verify that status is "0", not [None] as done below.  This peculiar
        # expected value is tested because we intercept the load() function,
        # because we lack a MongoDB emulator.
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], None)

    def test_bad_metadata_file(self):
        """Ingest Pipeline should not succeed for misformatted metadata file
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cell_metadata',
            '--cell-metadata-file',
            '../tests/data/metadata_bad.txt',
            '--study-accession',
            'SCP123',
            '--ingest-cell-metadata',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_bad_metadata_file_contains_coordinates(self):
        """Ingest Pipeline should not succeed for metadata file containing
        coordinates
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cell_metadata',
            '--cell-metadata-file',
            '../tests/data/metadata_bad_contains_coordinates.txt',
            '--study-accession',
            'SCP123',
            '--ingest-cell-metadata',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_good_cluster_file(self):
        """Ingest Pipeline should succeed for properly formatted cluster file
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cluster',
            '--cluster-file',
            '../tests/data/cluster_example.txt',
            '--ingest-cluster',
            '--name',
            'cluster1',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        # TODO:
        # After integrating MongoDB emulator (SCP-2000), refactor this test to
        # verify that status is "0", not [None] as done below.  This peculiar
        # expected value is tested because we intercept the load() function,
        # because we lack a MongoDB emulator.
        self.assertEqual(len(status), 1)
        self.assertEqual(status[0], None)

    def test_bad_cluster_file(self):
        """Ingest Pipeline should fail for misformatted cluster file
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cluster',
            '--cluster-file',
            '../tests/data/cluster_bad.txt',
            '--ingest-cluster',
            '--name',
            'cluster1',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    def test_bad_cluster_missing_coordinate_file(self):
        """Ingest Pipeline should fail for missing coordinate in cluster file
        """
        args = [
            '--study-id',
            '5d276a50421aa9117c982845',
            '--study-file-id',
            '5dd5ae25421aa910a723a337',
            'ingest_cluster',
            '--cluster-file',
            '../tests/data/cluster_bad_missing_coordinate.txt',
            '--ingest-cluster',
            '--name',
            'cluster1',
        ]
        ingest, arguments, status, status_cell_metadata = self.setup_ingest(args)

        with self.assertRaises(SystemExit) as cm:
            exit_pipeline(ingest, status, status_cell_metadata, arguments)
        self.assertEqual(cm.exception.code, 1)

    # def test_ingest_loom(self):
    #     """Ingest Pipeline should extract and transform loom files
    #     """
    #
    #     args = [
    #         '--study-id',
    #         '5d276a50421aa9117c982845',
    #         '--study-file-id',
    #         '5dd5ae25421aa910a723a337',
    #         'ingest_expression',
    #         '--taxon-name',
    #         'Homo Sapiens',
    #         '--taxon-common-name',
    #         'human',
    #         '--ncbi-taxid',
    #         '9606',
    #         '--genome-assembly-accession',
    #         'GCA_000001405.15',
    #         '--genome-annotation',
    #         'Ensemble 94',
    #         '--matrix-file',
    #         '../tests/data/test_loom.loom',
    #         '--matrix-file-type',
    #         'loom',
    #     ]
    #
    #     ingest = self.setup_ingest(args)
    #
    #     model = ingest.load_args[0]
    #
    #     # Verify that 25 gene models were passed into load method
    #     num_models = len(models)
    #     expected_num_models = 10
    #     self.assertEqual(num_models, expected_num_models)
    #
    #     # Verify that the first gene model looks as expected
    #     mock_dir = 'loom'
    #     model, expected_model = get_nth_gene_models(0, models, mock_dir)
    #     self.assertEqual(model, expected_model)


if __name__ == '__main__':
    unittest.main()
