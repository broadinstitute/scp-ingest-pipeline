import sys
import unittest
from argparse import ArgumentTypeError

sys.path.append("../ingest")
from cli_parser import is_valid_uuid


class TestAnnotations(unittest.TestCase):
    def test_is_valid_uuid(self):
        self.assertRaises(ArgumentTypeError, is_valid_uuid, "abc")
        self.assertFalse(is_valid_uuid(""))

        user_metrics_uuid = "123e4567-e89b-12d3-a456-426614174000"
        self.assertEqual(is_valid_uuid(user_metrics_uuid), user_metrics_uuid)
