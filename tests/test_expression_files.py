import sys
import unittest
from unittest.mock import MagicMock
from bson.objectid import ObjectId

sys.path.append("../ingest")
from expression_files.expression_files import GeneExpression


class TestExpressionFiles(unittest.TestCase):
    def test_has_unique_cells(self):
        client_mock = MagicMock()
        header = ["GENE", "foo", "foo2", "foo3"]

        client_mock["data_arrays"].find.return_value = [
            {"values": ["foo3", "foo4", "foo5"]},
            {"values": ["foo6", "foo7", "foo8"]},
        ]
        self.assertFalse(
            GeneExpression.has_unique_cells(header, ObjectId(), client_mock)
        )

        client_mock["data_arrays"].find.return_value = [
            {"values": ["foo4", "foo5", "foo6"]},
            {"values": ["foo9", "foo7", "foo8"]},
        ]
        self.assertTrue(
            GeneExpression.has_unique_cells(header, ObjectId(), client_mock)
        )

        client_mock["data_arrays"].find.return_value = []
        self.assertTrue(
            GeneExpression.has_unique_cells(header, ObjectId(), client_mock)
        )
