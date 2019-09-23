"""Tests for metadata validation

These tests verify that metadata files are checked against metadata convention,
ontology terms are validated against an external source, and tsv metadata files
conform to SCP metadata file format requirements.

PREREQUISITEs
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
and Firestore emulator must be running, see PR26 for instructions
(https://github.com/broadinstitute/scp-ingest-pipeline/pull/26)

Note: When CI environment moves to Python 3.7, tests may break due to minor
differences in how the reference issues are serialized

# Run all tests in a manner that shows report_issues output
python3 test_validate_metadata.py

"""

import sys
import unittest
import json

sys.path.append("../ingest")
sys.path.append("../ingest/validation")

from validate_metadata import (
    create_parser,
    report_issues,
    collect_jsonschema_errors,
    validate_schema,
    CellMetadata,
    validate_collected_ontology_data,
)


class TestValidateMetadata(unittest.TestCase):
    def setup_metadata(self, args):
        args_list = args.split(" ")
        args = create_parser().parse_args(args_list)
        with open(args.convention, "r") as f:
            convention = json.load(f)
        filetsv = args.input_metadata
        metadata = CellMetadata(filetsv, "1234abc", "SCP1")
        metadata.validate_format()
        return (metadata, convention)

    def teardown_metadata(self, metadata):
        metadata.file_handle.close()

    def test_header_format(self):
        """Header rows of metadata file should conform to standard
        """

        args = "../tests/data/AMC_v0.8.json " "../tests/data/error_headers.tsv"
        metadata = self.setup_metadata(args)[0]
        self.assertFalse(metadata.validate_header_keyword())
        self.assertIn(
            "Error: Metadata file header row malformed, missing NAME",
            metadata.issues['error']["format"].keys(),
            "Missing NAME keyword should fail format validation",
        )

        self.assertFalse(metadata.validate_type_keyword())
        self.assertIn(
            "Error: Metadata file TYPE row malformed, missing TYPE",
            metadata.issues['error']["format"].keys(),
            "Missing TYPE keyword should fail format validation",
        )

        self.assertFalse(
            metadata.validate_type_annotations(),
            "Invalid type annotations should fail format validation",
        )

        self.assertFalse(
            metadata.validate_unique_header(),
            "Duplicate headers should fail format validation",
        )

        self.assertFalse(
            metadata.validate_against_header_count(),
            "Mismatch between header and type annotation count "
            "should fail format validation",
        )

        self.assertTrue(
            report_issues(metadata), "Invalid metadata content should report issues"
        )

        self.teardown_metadata(metadata)

    def test_convention_content(self):
        """Metadata convention should be valid jsonschema
            """

        args = "../tests/data/AMC_invalid.json " "../tests/data/metadata_valid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.assertIsNone(
            validate_schema(convention, metadata),
            "Invalid metadata schema should be detected",
        )
        self.teardown_metadata(metadata)

    def test_valid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v0.8.json " "../tests/data/metadata_valid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid metadata content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_invalid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v1.1.00.json " "../tests/data/metadata_invalid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.maxDiff = None
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        self.assertTrue(
            report_issues(metadata), "Valid metadata content should not elicit error"
        )
        validate_collected_ontology_data(metadata, convention)
        # reference errors tests for:
        #   missing required property "sex"
        #   missing dependency for non-required property "ethinicity"
        #   missing value for non-required property "is_living"
        #   value provided not in enumerated list for "sample_type"
        #   value provided not a number for "organism_age"
        reference_file = open("../tests/data/metadata_invalid.json", "r")
        reference_issues = json.load(reference_file)
        reference_file.close()
        self.assertEqual(
            metadata.issues,
            reference_issues,
            "Metadata validation issues do not match reference issues",
        )

        self.teardown_metadata(metadata)

    def test_valid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v1.1.00.json " "../tests/data/ontology_valid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        validate_collected_ontology_data(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_invalid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v1.1.00.json " "../tests/data/ontology_invalid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.maxDiff = None
        metadata.validate_format()
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        validate_collected_ontology_data(metadata, convention)
        # reference errors tests for:
        #   invalid ontology shortname CELL for cell_type
        #   invalid ontologyID UBERON_1000331 for organ__ontology_label
        #   invalid ontology label "homo sapien" for species__ontology_label
        #     with species ontologyID of "NCBITaxon_9606"
        reference_file = open("../tests/data/ontology_invalid.json", "r")
        reference_issues = json.load(reference_file)
        reference_file.close()
        self.assertEqual(
            metadata.issues,
            reference_issues,
            "Ontology validation issues do not match reference issues",
        )

        self.teardown_metadata(metadata)


if __name__ == "__main__":
    unittest.main()
