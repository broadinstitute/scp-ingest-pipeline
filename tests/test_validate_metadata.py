import sys
import unittest
import json

sys.path.append("../ingest")
sys.path.append("../ingest/validation")

from validate_metadata import (
    create_parser,
    report_errors,
    process_metadata_content,
    validate_schema,
    CellMetadata,
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
            metadata.errors["format"],
            "Missing NAME keyword should fail format validation",
        )

        self.assertFalse(metadata.validate_type_keyword())
        self.assertIn(
            "Error:  Metadata file TYPE row malformed, missing TYPE",
            metadata.errors["format"],
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
            report_errors(metadata), "Invalid metadata content should report errors"
        )

        self.teardown_metadata(metadata)

    def test_convention_content(self):
        """Metadata convention should be valid jsonschema
            """

        args = "../tests/data/AMC_invalid.json " "../tests/data/metadata_valid.tsv"
        metadata, convention = self.setup_metadata(args)
        self.assertIsNone(
            validate_schema(convention), "Invalid metadata schema should be detected"
        )
        self.teardown_metadata(metadata)

    def test_valid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v0.8.json " "../tests/data/metadata_valid.tsv"
        metadata, convention = self.setup_metadata(args)
        metadata_valid = metadata.validate_format()
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        if metadata_valid:
            process_metadata_content(metadata, convention)
        self.assertFalse(
            report_errors(metadata), "Valid metadata content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_invalid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
            """
        args = "../tests/data/AMC_v0.8.json " "../tests/data/metadata_invalid.tsv"
        metadata, convention = self.setup_metadata(args)

        metadata_valid = metadata.validate_format()
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        if metadata_valid:
            process_metadata_content(metadata, convention)
        self.assertTrue(
            report_errors(metadata), "Valid metadata content should not elicit error"
        )

        # reference errors tests for:
        #   missing required property "sex"
        #   missing dependency for non-required property "ethinicity"
        #   missing value for non-required property "is_living"
        #   value provided not in enumerated list for "sample_type"
        #   value provided not a number for "organism_age"
        reference_file = open("../tests/data/metadata_invalid.json", "r")
        reference_errors = json.load(reference_file)
        reference_file.close()
        self.assertEqual(
            metadata.errors,
            reference_errors,
            "Metadata validation errors do not match reference errors",
        )
        self.assertEqual(
            metadata.errors["format"],
            reference_errors["format"],
            "Expected duplicate cellID error does not match reference",
        )

        self.teardown_metadata(metadata)


if __name__ == "__main__":
    unittest.main()
