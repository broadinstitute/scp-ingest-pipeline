"""Tests for metadata validation

These tests verify that metadata files are checked against metadata convention,
ontology terms are validated against an external source, and tsv metadata files
conform to SCP metadata file format requirements.

PREREQUISITES
Spin up Python 3.6 virtualenv, install Python dependencies in requirements.txt
and Firestore emulator must be running, see PR26 for instructions
(https://github.com/broadinstitute/scp-ingest-pipeline/pull/26)

# Run all tests in a manner that shows report_issues output
python3 test_validate_metadata.py

"""

import sys
import unittest
import json

sys.path.append('../ingest')
sys.path.append('../ingest/validation')

from validate_metadata import (
    create_parser,
    report_issues,
    collect_jsonschema_errors,
    validate_schema,
    CellMetadata,
    validate_collected_ontology_data,
    validate_input_metadata,
)


class TestValidateMetadata(unittest.TestCase):
    def setup_metadata(self, args):
        args_list = args.split(' ')
        args = create_parser().parse_args(args_list)
        with open(args.convention) as f:
            convention = json.load(f)
        filetsv = args.input_metadata
        metadata = CellMetadata(filetsv, '1234abc', 'SCP1')
        metadata.validate_format()
        print(f"Format is corrrect {metadata.validate_format()}")
        return (metadata, convention)

    def teardown_metadata(self, metadata):
        metadata.file_handle.close()

    def test_header_format(self):
        """Header rows of metadata file should conform to standard
        """

        args = '../tests/data/AMC_v1.1.2.json ../tests/data/error_headers_v1.1.1.tsv'
        metadata = self.setup_metadata(args)[0]
        self.assertFalse(metadata.validate_header_keyword())
        self.assertIn(
            'Malformed metadata file header row, missing NAME. (Case Sensitive)',
            metadata.issues['error']['format'].keys(),
            'Missing NAME keyword should fail format validation',
        )

        self.assertFalse(metadata.validate_type_keyword())
        self.assertIn(
            'Malformed metadata TYPE row, missing TYPE. (Case Sensitive)',
            metadata.issues['error']['format'].keys(),
            'Missing TYPE keyword should fail format validation',
        )

        self.assertFalse(
            metadata.validate_type_annotations(),
            'Invalid type annotations should fail format validation',
        )

        self.assertFalse(
            metadata.validate_unique_header(),
            'Duplicate headers should fail format validation',
        )

        self.assertTrue(
            report_issues(metadata), 'Invalid metadata content should report issues'
        )

        self.teardown_metadata(metadata)

    def test_convention_content(self):
        """Metadata convention should be valid jsonschema
        """

        args = '../tests/data/AMC_invalid.json ../tests/data/valid_v1.1.1.tsv'
        metadata, convention = self.setup_metadata(args)
        self.assertIsNone(
            validate_schema(convention, metadata),
            'Invalid metadata schema should be detected',
        )
        self.teardown_metadata(metadata)

    def test_valid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
        """
        # Note: this input metadata file does not have array-based metadata
        # is compatible with v1.1.2 but not v1.1.3 (missing sampleID and donorID)
        args = '../tests/data/AMC_v1.1.2.json ../tests/data/valid_v1.1.1.tsv'
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), 'Valid metadata headers should not elicit error'
        )
        collect_jsonschema_errors(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid metadata content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_invalid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
        """
        args = '../tests/data/AMC_v1.1.2.json ../tests/data/invalid_metadata_v1.1.1.tsv'
        metadata, convention = self.setup_metadata(args)
        self.maxDiff = None
        self.assertTrue(
            metadata.validate_format(), 'Valid metadata headers should not elicit error'
        )
        collect_jsonschema_errors(metadata, convention)
        self.assertTrue(
            report_issues(metadata), 'Valid metadata content should not elicit error'
        )
        validate_collected_ontology_data(metadata, convention)
        # reference errors tests for:
        #   missing required property 'sex'
        #   missing dependency for non-required property 'ethinicity'
        #   missing value for non-required property 'is_living'
        #   value provided not in enumerated list for 'sample_type'
        #   value provided not a number for 'organism_age'
        reference_file = open('../tests/data/issues_metadata_v1.1.1.json')
        reference_issues = json.load(reference_file)
        reference_file.close()
        self.assertEqual(
            metadata.issues,
            reference_issues,
            'Metadata validation issues do not match reference issues',
        )
        self.teardown_metadata(metadata)

    def test_valid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
        """
        # Note: this input metadata file does not have array-based metadata
        # is compatible with v1.1.2 but not v1.1.3 (missing sampleID and donorID)
        args = '../tests/data/AMC_v1.1.2.json ../tests/data/valid_v1.1.1.tsv'
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), 'Valid metadata headers should not elicit error'
        )
        validate_input_metadata(metadata, convention)
        self.assertFalse(
            report_issues(metadata), 'Valid ontology content should not elicit error'
        )
        self.teardown_metadata(metadata)

    def test_invalid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
        """
        # Note: this input metadata file does not have array-based metadata
        # is compatible with v1.1.2 but not v1.1.3 (missing sampleID and donorID)
        args = '../tests/data/AMC_v1.1.2.json ../tests/data/invalid_ontology_v1.1.1.tsv'
        metadata, convention = self.setup_metadata(args)
        self.maxDiff = None
        self.assertTrue(
            metadata.validate_format(), 'Valid metadata headers should not elicit error'
        )
        validate_input_metadata(metadata, convention)
        # reference errors tests for:
        #   empty cell for cell_type entry (convention and ontology errors)
        #   empty cell for entry of geographical_region and its label
        #       (convention and ontology errors)
        #   improper syntax (lack of _ or :) for EFO0008919
        #       (convention and ontology errors)
        #   invalid ontology shortname CELL for cell_type
        #   invalid ontology label 'homo sapien' for species__ontology_label
        #     with species ontologyID of 'NCBITaxon_9606'
        #   invalid ontologyID of 'NCBITaxon_9606' for geographical_region
        #   invalid ontologyID UBERON_1000331 for organ__ontology_label
        reference_file = open('../tests/data/issues_ontology_v1.1.1.json')
        reference_issues = json.load(reference_file)

        self.assertEqual(
            metadata.issues,
            reference_issues,
            'Ontology validation issues do not match reference issues',
        )
        reference_file.close()
        self.teardown_metadata(metadata)

    # def test_valid_array_content(self):
    #     """array-based metadata should conform to convention requirements
    #     """
    #     args = '../tests/data/AMC_v1.1.3.json ../tests/data/valid_array_v1.1.3.tsv'
    #     metadata, convention = self.setup_metadata(args)
    #     self.assertTrue(
    #         metadata.validate_format(), 'Valid metadata headers should not elicit error'
    #     )
    #     validate_input_metadata(metadata, convention)
    #     self.assertFalse(
    #         report_issues(metadata), 'Valid ontology content should not elicit error'
    #     )
    #     # valid array data emits one warning message for disease__time_since_onset__unit
    #     # because no ontology label supplied in metadata file for the unit ontology
    #     reference_file = open('../tests/data/issues_warn_v1.1.2.json')
    #     reference_issues = json.load(reference_file)
    #     reference_file.close()
    #     self.assertEqual(
    #         metadata.issues,
    #         reference_issues,
    #         'Metadata validation issues do not match reference issues',
    #     )
    #     self.teardown_metadata(metadata)
    #
    # def test_invalid_array_content(self):
    #     """array-based metadata should conform to convention requirements
    #     """
    #     args = '../tests/data/AMC_v1.1.3.json ../tests/data/invalid_array_v1.1.3.tsv'
    #     metadata, convention = self.setup_metadata(args)
    #     self.assertTrue(
    #         metadata.validate_format(), 'Valid metadata headers should not elicit error'
    #     )
    #     validate_input_metadata(metadata, convention)
    #     # reference errors tests for:
    #     # conflict between convention type and input metadata type annotation
    #     #     group instead of numeric: organism_age
    #     #     numeric instead of group: sample_type
    #     # invalid array-based metadata type: disease__time_since_onset
    #     # invalid boolean value: disease__treated
    #     # non-uniform unit values: organism_age__unit
    #     # missing ontology ID or label for non-required metadata: ethnicity
    #     reference_file = open('../tests/data/issues_array_v1.1.2.json')
    #     reference_issues = json.load(reference_file)
    #     reference_file.close()
    #     self.assertEqual(
    #         metadata.issues,
    #         reference_issues,
    #         'Metadata validation issues do not match reference issues',
    #     )
    #     self.teardown_metadata(metadata)


if __name__ == '__main__':
    unittest.main()
