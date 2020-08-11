import sys
import unittest
from unittest.mock import MagicMock
from bson.objectid import ObjectId

sys.path.append("../ingest")
from expression_files.expression_files import GeneExpression


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
