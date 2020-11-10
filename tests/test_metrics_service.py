import unittest
from unittest.mock import patch

import pytest
from unittest.mock import MagicMock
import json
from test_ingest import IngestTestCase
import config
from monitoring.metrics_service import MetricsService


# Mocked function that replaces MetricsService.post_event
def mock_post_event(props):
    expected_props = {
        "event": "ingest-pipeline:expression:ingest",
        "properties": {
            "studyAccession": "SCP123",
            "fileName": "File_name.txt",
            "fileType": "Expression matrix",
            "fileSize": 400,
            "appId": "single-cell-portal",
            "functionName": "ingest_expression",
            "distinct_id": MetricsService.user_id,
        },
    }
    props = json.loads(props)
    del props["properties"]["perfTime"]
    assert props == expected_props


class MetricsServiceTestCase(unittest.TestCase):
    client_values = {}
    client_values["study_accessions"] = MagicMock()
    client_values["study_accessions"].find.return_value = [{"accession": "SCP123"}]
    client_values["study_files"] = MagicMock()
    client_values["study_files"].find.return_value = [
        {
            "file_type": "Expression matrix",
            "upload_file_size": 400,
            "name": "File_name.txt",
        }
    ]
    # Build client mock with functions and return values needed to query
    client_mock = MagicMock()
    client_mock.__getitem__.side_effect = client_values.__getitem__

    @pytest.mark.usefixtures("delete_testing")
    @patch("expression_files.DenseIngestor.execute_ingest")
    @patch("config.MONGO_CONNECTION")
    @patch(
        "monitoring.metrics_service.MetricsService.post_event",
        side_effect=mock_post_event,
    )
    def test_ingest_dense_matrix(
        self, mock_execute_ingest, mock_MONGO_CONNECTION, mock_post_event
    ):
        mock_MONGO_CONNECTION._client = MetricsServiceTestCase.client_mock
        args = [
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
        # Initialize global variables
        config.init("5d276a50421aa9117c982845", "5dd5ae25421aa910a723a337")
        config.set_parent_event_name("ingest-pipeline:expression:ingest")
        IngestTestCase.execute_ingest(args)
        metrics_model = config.get_metric_properties()
        # Log Mixpanel events
        MetricsService.log(config.get_parent_event_name(), metrics_model)
