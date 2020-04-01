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
        project = args.project
        version = args.version
        return (project, version)

    def test_serialize_convention(self):
        """Verify known good input generates the expected metadata convention json
        """

        args = 'test 1.0.0'
        project, version = self.setup_metadata(args)
        self.maxDiff = None
        input_tsv = '../tests/data/test_convention.tsv'
        output_fullpath = 'output_convention.json'
        schema_info = build_schema_info(project, version)
        convention = serialize_convention(schema_info, input_tsv)
        write_schema(convention, output_fullpath)
        with open("output_convention.json", "r") as result:
            generated_convention = json.load(result)
        os.remove("output_convention.json")

        with open("../tests/data/test_convention.json", "r") as reference_file:
            reference_convention = json.load(reference_file)
        self.assertEqual(
            generated_convention,
            reference_convention,
            "Generated convention not match reference convention",
        )


if __name__ == "__main__":
    unittest.main()
