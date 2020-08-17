import sys
import unittest
from unittest.mock import MagicMock, patch, call
from bson.objectid import ObjectId


sys.path.append("../ingest")
from expression_files.expression_files import GeneExpression
from ingest_files import DataArray
from typing import Dict


class TestExpressionFiles(unittest.TestCase):
    def test_check_unique_cells(self):
        client_mock = MagicMock()
        header = ["GENE", "foo", "foo2", "foo3"]

        client_mock["data_arrays"].find.return_value = [
            {"values": ["foo3", "foo4", "foo5"]},
            {"values": ["foo6", "foo7", "foo8"]},
        ]
        with self.assertRaises(ValueError) as cm:
            GeneExpression.check_unique_cells(header, ObjectId(), client_mock)
        self.assertTrue("contains 1 cells" in str(cm.exception))
        self.assertTrue("foo3" in str(cm.exception))

        header = ["GENE", "foo", "foo3", "foo2", "foo6"]
        with self.assertRaises(ValueError) as cm:
            GeneExpression.check_unique_cells(header, ObjectId(), client_mock)
        self.assertTrue("contains 2 cells" in str(cm.exception))
        self.assertTrue("foo3" in str(cm.exception))
        self.assertTrue("foo6" in str(cm.exception))

        header = ["GENE", "foo", "foo2"]
        self.assertTrue(
            GeneExpression.check_unique_cells(header, ObjectId(), client_mock)
        )

        client_mock["data_arrays"].find.return_value = []
        self.assertTrue(
            GeneExpression.check_unique_cells(header, ObjectId(), client_mock)
        )

    def test_create_data_array(self):
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
            TypeError, [GeneExpression.create_data_array], "bad_position_arg", **kwargs
        )
        self.assertRaises(
            TypeError, [GeneExpression.create_data_array], hi="bad_kwarg", **kwargs
        )
        actual_da: Dict = next(GeneExpression.create_data_array(**kwargs))
        self.assertEqual(DataArray(**kwargs).__dict__, actual_da)

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

    @patch("expression_files.expression_files.GeneExpression.insert")
    def test_load(self, mock_insert):
        client_mock = MagicMock()
        study_id = ObjectId()
        study_file_id = ObjectId()
        expression_matrix = GeneExpression(
            "../tests/data/dense_matrix_19_genes_1000_cells.txt",
            study_file_id,
            study_id,
        )
        expression_matrix.mongo_connection = client_mock
        gene_docs = [
            {"name": "foo_study"},
            {"study_file_id": study_file_id},
            {"study_id": study_id},
        ]
        da_docs = [
            {"phi": study_file_id},
            {"fye": study_id},
            {"foo": study_id},
            {"phom": study_id},
        ]
        expected_calls = [
            # GeneExpression.insert() possible arguments
            call(gene_docs, expression_matrix.COLLECTION_NAME, client_mock),
            call(da_docs, "data_arrays", client_mock),
        ]

        expression_matrix.load(gene_docs, da_docs)
        expression_matrix.insert.assert_has_calls(expected_calls, any_order=False)

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
