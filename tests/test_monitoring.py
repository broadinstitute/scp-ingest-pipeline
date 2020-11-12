import unittest
from unittest.mock import patch

from test_ingest import IngestTestCase
from mock_data.client import mock_client
import config

import pytest


class MonitorTestCase(unittest.TestCase):
    ingest_args = [
        "--study-id",
        "5d276a50421aa9117c982845",
        "--study-file-id",
        "5dd5ae25421aa910a723a337",
        "ingest_expression",
        "--taxon-name",
        "Homo sapiens",
        "--taxon-common-name",
        "human",
        "--ncbi-taxid",
        "9606",
        "--genome-assembly-accession",
        "GCA_000001405.15",
        "--genome-annotation",
        "Ensembl 94",
        "--matrix-file",
        "gs://fake-bucket/tests/data/dense_matrix_19_genes_1000_cells.txt",
        "--matrix-file-type",
        "dense",
    ]

    @patch("expression_files.DenseIngestor.execute_ingest")
    @patch("config.MONGO_CONNECTION")
    @patch("monitoring.mixpanel_log.custom_metric")
    def test_guard_decorator_testing_not_set(
        self, mock_execute_ingest, mock_MONGO_CONNECTION, mock_custom_metric
    ):
        mock_MONGO_CONNECTION._client = mock_client
        # Initialize global variables
        config.init("5d276a50421aa9117c982845", "5dd5ae25421aa910a723a337")
        IngestTestCase.execute_ingest(MonitorTestCase.ingest_args)
        self.assertTrue(mock_custom_metric.called)

    @pytest.mark.usefixtures("delete_testing")
    @patch("expression_files.DenseIngestor.execute_ingest")
    @patch("config.MONGO_CONNECTION")
    @patch("monitoring.mixpanel_log.custom_metric")
    def test_guard_decorator_testing_set(
        self, mock_execute_ingest, mock_MONGO_CONNECTION, mock_custom_metric
    ):
        mock_MONGO_CONNECTION._client = MonitorTestCase.client_mock
        # Initialize global variables
        config.init("5d276a50421aa9117c982845", "5dd5ae25421aa910a723a337")
        IngestTestCase.execute_ingest(MonitorTestCase.ingest_args)
        self.assertFalse(mock_custom_metric.called)
