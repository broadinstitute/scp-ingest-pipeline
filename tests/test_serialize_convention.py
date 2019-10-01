"""Tests for serialize_convention

This test verifies that a known good input metadata convention TSV file
generates the expected metadata convention output

PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt

# Run test
python3 test_serialize_convention.py

"""

import sys
import unittest
import json
import os

sys.path.append("../scripts")

from serialize_convention import (
    create_parser,
    build_schema_info,
    serialize_convention,
    write_schema,
)


class TestValidateMetadata(unittest.TestCase):
    def setup_metadata(self, args):
        args_list = args.split(" ")
        args = create_parser().parse_args(args_list)
        input_convention = args.input_convention
        label = args.label
        project = args.project
        return (input_convention, label, project)

    def test_serialize_convention(self):
        """Verify known good input generates the expected metadata convention json
        """

        args = "test ../tests/data/test_convention.tsv"
        input_convention, label, project = self.setup_metadata(args)
        self.maxDiff = None
        schema_info = build_schema_info(project)
        convention = serialize_convention(schema_info, input_convention)
        write_schema(convention, input_convention, label, ".")
        with open("test_convention.json", "r") as result:
            generated_convention = json.load(result)
        os.remove("test_convention.json")

        with open("../tests/data/test_convention.json", "r") as reference_file:
            reference_convention = json.load(reference_file)
        self.assertEqual(
            generated_convention,
            reference_convention,
            "Generated convention not match reference convention",
        )


if __name__ == "__main__":
    unittest.main()
