import sys
import unittest
from unittest.mock import patch, MagicMock
from bson.objectid import ObjectId

from pymongo.errors import AutoReconnect, BulkWriteError

# Dense models
from mock_data.expression.dense_matrices.nineteen_genes_100k_cell_models import (
    nineteen_genes_100k_cell_models,
)

# MTX models
from mock_data.expression.matrix_mtx.AB_toy_data_toy_models import (
    AB_toy_data_toy_data_models,
)
from mock_data.expression.matrix_mtx.one_column_feature_file_data_models import (
    one_column_feature_file_data_models,
)
from mock_data.expression.matrix_mtx.AB_toy_data_toy_unsorted_mtx_models import (
    AB_toy_data_toy_unsorted_mtx_models,
)

sys.path.append("../ingest")
from expression_files.expression_files import GeneExpression
from ingest_files import DataArray
from typing import Dict


def mock_load_genes_batched(documents, collection_name):
    """Enables overwriting of GeneExpression.load() with this placeholder.
    GeneExpression.load() is called multiple times. This method will verify
    data array documents are batched correctly everytime a call is made to
    GeneExpression.load()
    """
    # Assures transform function batches data array creation correctly
    if collection_name == DataArray.COLLECTION_NAME:
        assert GeneExpression.DATA_ARRAY_BATCH_SIZE >= len(documents)


def mock_expression_load(self, *args):
    """Enables overwriting of GeneExpression.load() with this placeholder.
    GeneExpression.load() is called multiple times. This method will verify
    models in the arguments have the expected values.
    """
    documents = args[0]
    collection_name = args[1]
    if not self.test_models:
        model_name = documents[0]["name"]
        if collection_name == GeneExpression.COLLECTION_NAME:
            if model_name in AB_toy_data_toy_data_models["gene_models"]:
                self.test_model = AB_toy_data_toy_data_models
            if model_name in AB_toy_data_toy_unsorted_mtx_models["gene_models"]:
                self.test_models = AB_toy_data_toy_unsorted_mtx_models
            if model_name in one_column_feature_file_data_models["gene_models"]:
                self.test_models = one_column_feature_file_data_models
            if model_name in nineteen_genes_100k_cell_models["gene_models"][model_name]:
                self.test_models = nineteen_genes_100k_cell_models
        if collection_name == DataArray.COLLECTION_NAME:
            if model_name in AB_toy_data_toy_data_models["data_arrays"]:
                self.test_models = AB_toy_data_toy_data_models
            if model_name in AB_toy_data_toy_unsorted_mtx_models["data_arrays"]:
                self.test_models = AB_toy_data_toy_unsorted_mtx_models
            if model_name in one_column_feature_file_data_models["data_arrays"]:
                self.test_models = one_column_feature_file_data_models
            if model_name in nineteen_genes_100k_cell_models["data_arrays"]:
                self.test_models = nineteen_genes_100k_cell_models
    # _id and linear_data_id are unique identifiers and can not be predicted
    # so we exclude it from the comparison
    for document in documents:
        model_name = document["name"]
        if collection_name == GeneExpression.COLLECTION_NAME:
            expected_model = self.test_models["gene_models"][model_name]
            del document["_id"]
            # Sometimes linear_data_id is not deleted in mock models in cases where mocked models are large
            if "linear_data_id" in expected_model:
                del expected_model["linear_data_id"]
            assert document == expected_model
        if collection_name == DataArray.COLLECTION_NAME:
            expected_model = self.test_models["data_arrays"][model_name]
            del document["linear_data_id"]
            # Sometimes linear_data_id is not deleted in mock models in cases where mocked models are large
            if "linear_data_id" in expected_model:
                del expected_model["linear_data_id"]
            assert document == expected_model
    self.models_processed += len(documents)


client_values = {}
client_values["data_arrays"] = MagicMock()
client_values["data_arrays"].find.return_value = [
    {"values": ["foo3", "foo4", "foo5"]},
    {"values": ["foo6", "foo7", "foo8"]},
]
client_values["study_files"] = MagicMock()
client_values["study_files"].find.return_value = [
    {"_id": ObjectId("5f70abd6771a5b0de0cea0f0")},
    {"_id": ObjectId("5f70c1b2771a5b0de0cea0ed")},
]


class TestExpressionFiles(unittest.TestCase):
    STUDY_ID = ObjectId("5d276a50421aa9117c982845")
    STUDY_FILE_ID = ObjectId("5f70abd6771a5b0de0cea0f0")

    # Build client mock with functions and return values needed to query
    # Taken from https://stackoverflow.com/questions/10043965
    client_mock: MagicMock = MagicMock()
    client_mock.keys.return_value.__iter__.return_value = client_values.keys()
    client_mock.__getitem__.side_effect = lambda k: client_values[k]

    def setUp(self) -> None:
        # Reinitialize client mock before each test
        TestExpressionFiles.client_mock.keys.return_value.__iter__.return_value = (
            client_values.keys()
        )
        TestExpressionFiles.client_mock.__getitem__.side_effect = lambda k: client_values[
            k
        ]

    def test_check_unique_cells(self):
        with patch(
            "expression_files.expression_files.GeneExpression.get_cell_names_from_study_file_id",
            return_value=["foo3", "foo4", "foo5", "foo6", "foo7", "foo8"],
        ), self.assertRaises(ValueError) as cm:
            header = ["GENE", "foo", "foo2", "foo3"]
            GeneExpression.check_unique_cells(
                header, ObjectId(), ObjectId(), TestExpressionFiles.client_mock
            )
            self.assertTrue("contains 1 cells" in str(cm.exception))
            self.assertTrue("foo3" in str(cm.exception))

        # Duplicate cells but duplicate values don't appear in same 'value' array
        with patch(
            "expression_files.expression_files.GeneExpression.get_cell_names_from_study_file_id",
            return_value=["foo3", "foo4", "foo5", "foo6", "foo7", "foo8"],
        ), self.assertRaises(ValueError) as cm:
            header = ["GENE", "foo", "foo3", "foo2", "foo6"]
            GeneExpression.check_unique_cells(
                header, ObjectId(), ObjectId(), TestExpressionFiles.client_mock
            )
            self.assertTrue("contains 2 cells" in str(cm.exception))
            self.assertTrue("foo3" in str(cm.exception))
            self.assertTrue("foo6" in str(cm.exception))

        with patch(
            "expression_files.expression_files.GeneExpression.get_cell_names_from_study_file_id",
            return_value=None,
        ):
            # Cells are unique
            header = ["GENE", "foo", "foo2"]
            self.assertTrue(
                GeneExpression.check_unique_cells(
                    header,
                    TestExpressionFiles.STUDY_ID,
                    TestExpressionFiles.STUDY_FILE_ID,
                    TestExpressionFiles.client_mock,
                )
            )
        with patch(
            "expression_files.expression_files.GeneExpression.get_cell_names_from_study_file_id",
            return_value=["foo3", "foo4", "foo5", "foo6", "foo7", "foo8"],
        ):
            # Cells are unique
            header = ["GENE", "foo", "foo2"]
            self.assertTrue(
                GeneExpression.check_unique_cells(
                    header,
                    TestExpressionFiles.STUDY_ID,
                    TestExpressionFiles.STUDY_FILE_ID,
                    TestExpressionFiles.client_mock,
                )
            )

    @patch("expression_files.expression_files.GeneExpression.query_cells")
    def test_get_cell_names_from_study_file_id(self, mock_query_cells):
        with patch(
            "expression_files.expression_files.GeneExpression.get_study_expression_file_ids",
            return_value=[
                {"_id": ObjectId(), "study_file_id": ObjectId()},
                {"_id": ObjectId(), "study_file_id": ObjectId()},
            ],
        ):
            mock_query_cells.return_value = [
                {"values": ["foo3", "foo4", "foo5"]},
                {"values": ["foo6", "foo7", "foo8", "foo6", "foo7", "foo8"]},
            ]
            expected = [
                "foo3",
                "foo4",
                "foo5",
                "foo6",
                "foo7",
                "foo8",
                "foo6",
                "foo7",
                "foo8",
            ]
            cells = GeneExpression.get_cell_names_from_study_file_id(
                ObjectId(), ObjectId(), MagicMock()
            )
            self.assertEqual(cells, expected)

            mock_query_cells.return_value = []
            cells = GeneExpression.get_cell_names_from_study_file_id(
                ObjectId(), ObjectId(), MagicMock()
            )
            self.assertEqual(cells, None)

    def test_is_raw_count(self):
        client = MagicMock()
        client["study_files"].find.return_value = [
            {"expression_file_info": {"is_raw_counts": False}}
        ]
        self.assertFalse(
            GeneExpression.is_raw_count(
                TestExpressionFiles.STUDY_ID, TestExpressionFiles.STUDY_FILE_ID, client
            )
        )
        client["study_files"].find.return_value = [{}]
        self.assertFalse(
            GeneExpression.is_raw_count(
                TestExpressionFiles.STUDY_ID,
                TestExpressionFiles.STUDY_FILE_ID,
                TestExpressionFiles.client_mock,
            )
        )

        client["study_files"].find.return_value = [
            {"expression_file_info": {"is_raw_counts": True}}
        ]
        self.assertTrue(
            GeneExpression.is_raw_count(
                TestExpressionFiles.STUDY_ID, TestExpressionFiles.STUDY_FILE_ID, client
            )
        )

        client["study_files"].find.return_value = [
            {"expression_file_info": {"is_raw_counts": False}}
        ]
        self.assertFalse(
            GeneExpression.is_raw_count(
                TestExpressionFiles.STUDY_ID, TestExpressionFiles.STUDY_FILE_ID, client
            )
        )

    def test_get_study_expression_file_ids(self):
        RAW_COUNTS_QUERY = {
            "$and": [{"study_id": ObjectId("5d276a50421aa9117c982845")}],
            "file_type": {"$in": ["Expression Matrix", "MM Coordinate Matrix"]},
            "$nor": [{"_id": ObjectId("5f70abd6771a5b0de0cea0f0")}],
        }
        FIELD_NAME = {"_id": 1}
        COLLECTION_NAME = "study_files"

        # Study has study files with document expression_file_info

        with patch(
            "expression_files.expression_files.GeneExpression.is_raw_count",
            return_value=True,
        ):
            # Add query filter for is_raw_counts
            RAW_COUNTS_QUERY["$and"].append(
                {"expression_file_info.is_raw_counts": True}
            )
            GeneExpression.get_study_expression_file_ids(
                TestExpressionFiles.STUDY_ID,
                TestExpressionFiles.STUDY_FILE_ID,
                TestExpressionFiles.client_mock,
            )
            TestExpressionFiles.client_mock[COLLECTION_NAME].find.assert_called_with(
                RAW_COUNTS_QUERY, FIELD_NAME
            )

        # Return query to original state
        del RAW_COUNTS_QUERY["$and"][1]

        # When study file is not raw count
        with patch(
            "expression_files.expression_files.GeneExpression.is_raw_count",
            return_value=False,
        ):
            GeneExpression.get_study_expression_file_ids(
                TestExpressionFiles.STUDY_ID,
                TestExpressionFiles.STUDY_FILE_ID,
                TestExpressionFiles.client_mock,
            )
            TestExpressionFiles.client_mock[COLLECTION_NAME].find.assert_called_with(
                RAW_COUNTS_QUERY, FIELD_NAME
            )

    def test_create_data_arrays(self):
        _id = ObjectId()
        study_id = ObjectId()
        study_file_id = ObjectId()
        kwargs = {
            "name": "foo",
            "cluster_name": "foo_name",
            "array_type": "foo_type",
            "values": [1, 2, 3, 4, 5],
            "linear_data_type": "data_dtype",
            "linear_data_id": _id,
            "study_id": study_id,
            "study_file_id": study_file_id,
        }

        self.assertRaises(
            TypeError, [GeneExpression.create_data_arrays], "bad_position_arg", **kwargs
        )
        self.assertRaises(
            TypeError, [GeneExpression.create_data_arrays], hi="bad_kwarg", **kwargs
        )
        actual_data_array: Dict = next(GeneExpression.create_data_arrays(**kwargs))
        self.assertEqual(DataArray(**kwargs).__dict__, actual_data_array)

    def test_create_gene_model(self):
        _id = ObjectId()
        study_id = ObjectId()
        study_file_id = ObjectId()
        kwargs = {
            "name": "foo",
            "study_file_id": study_file_id,
            "study_id": study_id,
            "_id": _id,
        }
        self.assertRaises(
            TypeError, [GeneExpression.create_gene_model], "bad_position_arg", **kwargs
        )
        self.assertRaises(
            TypeError, [GeneExpression.create_gene_model], hi="bad_kwarg", **kwargs
        )
        actual_gene_model = GeneExpression.create_gene_model(**kwargs)
        expected_gene_model = GeneExpression.Model(
            searchable_name="foo", gene_id=None, **kwargs
        )
        self.assertEqual(actual_gene_model, expected_gene_model)

    def test_insert(self):

        client_mock = MagicMock()

        docs = [
            {"values": ["foo3", "foo4", "foo5"]},
            {"values": ["foo6", "foo7", "foo8"]},
        ]
        GeneExpression.insert(docs, "collection", client_mock)
        client_mock["collection"].insert_many.assert_called_with(docs, ordered=False)

        client_mock["collection"].insert_many.side_effect = ValueError("Foo")
        self.assertRaises(
            Exception, GeneExpression.insert, docs, "collection", client_mock
        )
        client_mock.reset_mock()

        # Test exponential back off for auto reconnect
        client_mock["collection"].insert_many.side_effect = AutoReconnect
        self.assertRaises(
            AutoReconnect, GeneExpression.insert, docs, "collection", client_mock
        )
        self.assertEqual(client_mock["collection"].insert_many.call_count, 5)
        client_mock.reset_mock()

        def raiseError(*args, **kwargs):
            raise BulkWriteError({"details": "foo"})

        # Test exponential back off for BulkWriteError
        client_mock["collection"].insert_many.side_effect = raiseError
        self.assertRaises(
            BulkWriteError, GeneExpression.insert, docs, "collection", client_mock
        )
        self.assertEqual(client_mock["collection"].insert_many.call_count, 5)
