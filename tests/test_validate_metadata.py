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
import os
from bson.objectid import ObjectId
import requests
from unittest.mock import patch
import io
import numpy as np
import pandas as pd
import pprint

sys.path.append("../ingest")
sys.path.append("../ingest/validation")

from cell_metadata import CellMetadata
from validate_metadata import (
    report_issues,
    collect_jsonschema_errors,
    validate_schema,
    validate_collected_ontology_data,
    collect_cell_for_ontology,
    validate_input_metadata,
    request_json_with_backoff,
    MAX_HTTP_ATTEMPTS,
    is_empty_string,
    is_label_or_synonym,
    replace_single_value_array,
    replace_synonym_in_multivalue_array,
)
from metadata_validation import create_parser


# do not attempt a request, but instead throw a request exception
def mocked_requests_get(*args, **kwargs):
    raise requests.exceptions.RequestException


class TestValidateMetadata(unittest.TestCase):
    def setup_metadata(self, args):
        args_list = args.split(" ")
        args = create_parser().parse_args(args_list)
        with open(args.convention) as f:
            convention = json.load(f)
        filetsv = args.input_metadata
        # set ObjectIDs to be recognizably artificial
        artificial_study_file_id = "addedfeed000000000000000"
        metadata = CellMetadata(
            filetsv,
            ObjectId("dec0dedfeed1111111111111"),
            ObjectId(artificial_study_file_id),
            "SCPtest",
            study_accession="SCPtest",
        )
        metadata.preprocess(is_metadata_convention=True)
        # The following line should become metadata.validate() as part of bugfix SCP-2756
        metadata.validate_format()
        print(f"Format is correct {metadata.validate_format()}")
        return (metadata, convention)

    def teardown_metadata(self, metadata):
        metadata.file_handle.close()
        try:
            os.remove("scp_validation_errors.txt")
            os.remove("scp_validation_warnings.txt")
            os.remove("errors.txt")
            os.remove("info.txt")
        except OSError:
            print("no file to remove")

    def test_comma_delimited_array(self):
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/has_commas_in_arrays.csv"
        )
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        validate_collected_ontology_data(metadata, convention)
        # reference errors tests for:
        # incorrectly delimited value in numeric column
        self.assertIn(
            "disease__time_since_onset: '12,2' in '12,2' does not match expected 'number' type.",
            metadata.issues["error"]["content"].keys(),
            "comma-delimited values for array metadata should fail",
        )
        # incorrectly delimited value in boolean column
        self.assertIn(
            "disease__treated: 'True,False' in 'True,False' does not match expected 'boolean' type.",
            metadata.issues["error"]["content"].keys(),
            "comma-delimited values for array metadata should fail",
        )
        # partial correctly delimited value in numeric column should still fail
        self.assertIn(
            "disease__time_since_onset: '3,1' in '36|3,1' does not match expected 'number' type.",
            metadata.issues["error"]["content"].keys(),
            "incorrectly delimited values for array metadata should fail",
        )

        self.teardown_metadata(metadata)

    def test_convention_content(self):
        """Metadata convention should be valid jsonschema
        """

        args = "--convention ../tests/data/AMC_invalid.json ../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt"
        metadata, convention = self.setup_metadata(args)
        self.assertIsNone(
            validate_schema(convention, metadata),
            "Invalid metadata schema should be detected",
        )
        self.teardown_metadata(metadata)

    def test_auto_filling_missing_labels(self):
        # note that the filename provided here is irrelevant -- we will be specifying row data ourselves
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt"
        )
        metadata, convention = self.setup_metadata(args)

        # handle empty string ontology label for required array metadata
        input_row = {
            "CellID": "test1",
            "disease": ["MONDO_0005015"],
            "disease__ontology_label": "",
        }
        expected_row = {
            "CellID": "test1",
            "disease": ["MONDO_0005015"],
            "disease__ontology_label": [],
        }
        updated_row = collect_cell_for_ontology(
            "disease", input_row, metadata, convention, True, True
        )
        self.assertEqual(
            expected_row,
            updated_row,
            "Row should not be altered if label for required ontology is missing",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {
                'disease: required column "disease__ontology_label" missing data': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handle missing ontology label column for required array metadata
        metadata, convention = self.setup_metadata(args)
        row = {"CellID": "test1", "disease": ["MONDO_0005015"]}
        updated_row = collect_cell_for_ontology(
            "disease", row, metadata, convention, True, True
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "disease": ["MONDO_0005015"],
                "disease__ontology_label": [],
            },
            updated_row,
            "Row should have column ontology_label added with value of empty array",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {
                'disease: required column "disease__ontology_label" missing data': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handles nan ontology label for required array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "disease": ["MONDO_0005015"],
            "disease__ontology_label": np.nan,
        }
        updated_row = collect_cell_for_ontology(
            "disease", row, metadata, convention, True, True
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "disease": ["MONDO_0005015"],
                "disease__ontology_label": [],
            },
            updated_row,
            "nan should be converted to empty array",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {
                'disease: required column "disease__ontology_label" missing data': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handle empty string ontology label for required non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "organ": "UBERON_0001913",
            "organ__ontology_label": "",
        }
        updated_row = collect_cell_for_ontology(
            "organ", row, metadata, convention, False, True
        )
        self.assertEqual(
            row,
            updated_row,
            "Row should not be altered if label for required ontology is missing",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {'organ: required column "organ__ontology_label" missing data': ["test1"]},
            "unexpected error reporting",
        )

        # handle missing ontology label column for required non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {"CellID": "test1", "organ": "UBERON_0001913"}
        updated_row = collect_cell_for_ontology(
            "organ", row, metadata, convention, False, True
        )
        self.assertEqual(
            {"CellID": "test1", "organ": "UBERON_0001913", "organ__ontology_label": ""},
            updated_row,
            "Row should have column ontology_label added with value of empty string",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {'organ: required column "organ__ontology_label" missing data': ["test1"]},
            "unexpected error reporting",
        )

        # handles nan ontology label for required non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "organ": "UBERON_0001913",
            "organ__ontology_label": np.nan,
        }
        updated_row = collect_cell_for_ontology(
            "organ", row, metadata, convention, False, True
        )
        self.assertEqual(
            {"CellID": "test1", "organ": "UBERON_0001913", "organ__ontology_label": ""},
            updated_row,
            "nan should be converted to empty string",
        )
        self.assertEqual(
            metadata.issues["error"]["content"],
            {'organ: required column "organ__ontology_label" missing data': ["test1"]},
            "unexpected error reporting",
        )

        # handle empty string ontology label for optional array metadata
        row = {
            "CellID": "test1",
            "ethnicity": ["HANCESTRO_0005"],
            "ethnicity__ontology_label": "",
        }
        updated_row = collect_cell_for_ontology(
            "ethnicity", row, metadata, convention, True, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "ethnicity": ["HANCESTRO_0005"],
                "ethnicity__ontology_label": ["European"],
            },
            updated_row,
            "Row should be updated to inject missing ontology label as array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'ethnicity: missing ontology label "HANCESTRO_0005" - using "European" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handle missing ontology label column for optional array metadata
        metadata, convention = self.setup_metadata(args)
        row = {"CellID": "test1", "ethnicity": ["HANCESTRO_0005"]}
        updated_row = collect_cell_for_ontology(
            "ethnicity", row, metadata, convention, True, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "ethnicity": ["HANCESTRO_0005"],
                "ethnicity__ontology_label": ["European"],
            },
            updated_row,
            "Row should be updated to inject missing ontology label as array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'ethnicity: missing ontology label "HANCESTRO_0005" - using "European" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handles nan ontology label for optional array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "ethnicity": ["HANCESTRO_0005"],
            "ethnicity__ontology_label": np.nan,
        }
        updated_row = collect_cell_for_ontology(
            "ethnicity", row, metadata, convention, True, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "ethnicity": ["HANCESTRO_0005"],
                "ethnicity__ontology_label": ["European"],
            },
            updated_row,
            "Row should be updated to inject missing ontology label as array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'ethnicity: missing ontology label "HANCESTRO_0005" - using "European" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handle empty string ontology label for optional non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "cell_type": "CL_0000066",
            "cell_type__ontology_label": "",
        }
        updated_row = collect_cell_for_ontology(
            "cell_type", row, metadata, convention, False, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "cell_type": "CL_0000066",
                "cell_type__ontology_label": "epithelial cell",
            },
            updated_row,
            "Row should be updated to inject missing ontology label as non-array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'cell_type: missing ontology label "CL_0000066" - using "epithelial cell" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handle missing ontology label column for optional non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {"CellID": "test1", "cell_type": "CL_0000066"}
        updated_row = collect_cell_for_ontology(
            "cell_type", row, metadata, convention, False, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "cell_type": "CL_0000066",
                "cell_type__ontology_label": "epithelial cell",
            },
            updated_row,
            "Row should be updated to inject missing ontology label as non-array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'cell_type: missing ontology label "CL_0000066" - using "epithelial cell" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handles nan ontology label for optional non-array metadata
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "cell_type": "CL_0000066",
            "cell_type__ontology_label": np.nan,
        }
        updated_row = collect_cell_for_ontology(
            "cell_type", row, metadata, convention, False, False
        )
        self.assertEqual(
            {
                "CellID": "test1",
                "cell_type": "CL_0000066",
                "cell_type__ontology_label": "epithelial cell",
            },
            updated_row,
            "Row should be updated to inject missing ontology label as non-array",
        )
        self.assertEqual(
            metadata.issues["warn"]["ontology"],
            {
                'cell_type: missing ontology label "CL_0000066" - using "epithelial cell" per EBI OLS lookup': [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

        # handles mismatch in item #s for optional array metadata and its label
        metadata, convention = self.setup_metadata(args)
        row = {
            "CellID": "test1",
            "ethnicity": ["HANCESTRO_0005", "HANCESTRO_0462"],
            "ethnicity__ontology_label": ["British"],
        }
        updated_row = collect_cell_for_ontology(
            "ethnicity", row, metadata, convention, True, False
        )
        self.assertEqual(
            row,
            updated_row,
            "Row should not be altered if mismatch in item #s for between array metadata and its label",
        )
        self.assertEqual(
            metadata.issues["error"]["ontology"],
            {
                "ethnicity: mismatched # of ethnicity and ethnicity__ontology_label values.": [
                    "test1"
                ]
            },
            "unexpected error reporting",
        )

    def test_valid_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
        """
        # def set
        # Note: this input metadata file does not have array-based metadata
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt"
        )
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        collect_jsonschema_errors(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid metadata content should not elicit error"
        )
        self.teardown_metadata(metadata)

        # invalid non-ontology content
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/invalid_metadata_v2.0.0.tsv"
        )
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
        #   duplicate CellIDs are not permitted
        self.assertIn(
            "Duplicate CellID(s) in metadata file.",
            metadata.issues["error"]["content"].keys(),
            "Duplicate CellID(s) in metadata file should result in error",
        )

        #   missing required property 'sex'
        self.assertIn(
            "'sex' is a required property",
            metadata.issues["error"]["convention"].keys(),
            "absence of required metadata should result in error",
        )
        #  ontology ID required if ontology_label provided even for non-required metadata
        self.assertIn(
            "'ethnicity' is a dependency of 'ethnicity__ontology_label'",
            metadata.issues["error"]["convention"].keys(),
            "ontology_label without ontology ID data results in error, even for non-required metadata",
        )
        #   expect time unit if organism age is provided
        self.assertIn(
            "'organism_age' is a dependency of 'organism_age__unit'",
            metadata.issues["error"]["convention"].keys(),
            "absence of units for age metadata should result in error",
        )
        #   value provided not a number for 'organism_age'
        self.assertIn(
            "organism_age: \"foo\" does not match expected type.",
            metadata.issues["error"]["content"].keys(),
            "ontology_label without ontology ID data results in error, even for non-required metadata",
        )
        self.teardown_metadata(metadata)

    def test_valid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
        """
        # Note: this input metadata file does not have array-based metadata
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt"
        )
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        validate_input_metadata(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_valid_multiple_ontologies_content(self):
        """Ontology metadata should conform to convention requirements
           Specifically tests that a term can be found in one of two accepted ontologies (e.g. disease in MONDO or PATO)
        """
        # Note: this input metadata file does not have array-based metadata
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_no_array_v2.0.0.txt"
        )
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        validate_input_metadata(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_invalid_ontology_content(self):
        """Ontology metadata should conform to convention requirements
        """
        # Note: this input metadata file does not have array-based metadata
        args = "--convention ../schema/alexandria_convention/alexandria_convention_schema.json ../tests/data/invalid_ontology_v2.0.0.tsv"
        metadata, convention = self.setup_metadata(args)
        self.maxDiff = None
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        validate_input_metadata(metadata, convention)
        # reference errors tests for:
        #   empty cell for cell_type entry that
        #   has cell_type__ontology_label value (convention error)
        self.assertIn(
            "'cell_type' is a dependency of 'cell_type__ontology_label'",
            metadata.issues["error"]["convention"].keys(),
        )
        #   improper syntax (lack of _ or :) for EFO0008919
        #       (convention and ontology errors)
        self.assertIn(
            "library_preparation_protocol: 'EFO0008919' does not match '^[-A-Za-z0-9]+[_:][-A-Za-z0-9]+'",
            metadata.issues["error"]["convention"].keys(),
        )
        self.assertIn(
            'library_preparation_protocol: Could not parse provided ontology id, "EFO0008919".',
            metadata.issues["error"]["ontology"].keys(),
        )
        #   invalid ontology shortname CELL for cell_type
        self.assertIn(
            "cell_type: No match found in EBI OLS for provided ontology ID: CELL_0000066",
            metadata.issues["error"]["ontology"].keys(),
        )
        #   invalid ontology label 'homo sapien' for species__ontology_label
        #     with species ontologyID of 'NCBITaxon_9606'
        self.assertIn(
            'species: input ontology_label "homo sapien" does not match EBI OLS lookup "Homo sapiens" for ontology id "NCBITaxon_9606".',
            metadata.issues["error"]["ontology"].keys(),
        )
        #   invalid ontologyID 'NCBITaxon_9606' for geographical_region
        self.assertIn(
            "geographical_region: No match found in EBI OLS for provided ontology ID: NCBITaxon_9606",
            metadata.issues["error"]["ontology"].keys(),
        )
        #   invalid ontologyID 'UBERON_1000331' for organ__ontology_label
        self.assertIn(
            "organ: No match found in EBI OLS for provided ontology ID: UBERON_1001913",
            metadata.issues["error"]["ontology"].keys(),
        )
        #   improper ontology IDs 'NCIT_C-43862' and 'NC-IT_C126538' for race__ontology_label
        self.assertIn(
            "race: No match found in EBI OLS for provided ontology ID: NCIT_C-43862",
            metadata.issues["error"]["ontology"].keys(),
        )
        self.assertIn(
            "race: No match found in EBI OLS for provided ontology ID: NC-IT_C126538",
            metadata.issues["error"]["ontology"].keys(),
        )
        self.teardown_metadata(metadata)

    def test_content(self):
        """Array-based metadata should conform to convention requirements
        """

        def set_up_test(test_file_name, bq=None):
            test_file_path = "data/annotation/metadata/convention/" + test_file_name
            args = "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            metadata, convention = self.setup_metadata(args + test_file_path)
            if bq:
                validate_input_metadata(metadata, convention, bq_json=True)
            else:
                validate_input_metadata(metadata, convention)

            return metadata

        # metadata validation stores warning messages which are available thru CLI
        # but not currently included in error email - too noisy
        # because no ontology label supplied in metadata file for the unit ontology
        metadata = set_up_test("valid_array_v2.1.2.txt")
        self.assertIn(
            'disease__time_since_onset__unit: missing ontology label "UO_0000035" - using "month" per EBI OLS lookup',
            metadata.issues["warn"]["ontology"].keys(),
        )
        self.assertIn(
            'organ_region: missing ontology label "MBA_000000714" - using "Orbital area" per Mouse Brain Atlas ontology',
            metadata.issues["warn"]["ontology"].keys(),
        )
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

        # cell_type__custom entries no longer required to have corresponding cell_type metadata
        metadata = set_up_test("valid_cell_type__custom_v2.2.1.txt", bq=True)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

        # Negative test cases
        metadata = set_up_test("invalid_array_v2.1.2.txt")
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        # reference errors tests for:
        # conflict between convention type and input metadata type annotation
        #     group instead of numeric: organism_age
        self.assertIn(
            'organism_age: \"group\" annotation in metadata file conflicts with metadata convention. Convention expects \"numeric\" values.',
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        #     numeric instead of group: biosample_type
        self.assertIn(
            'biosample_type: \"numeric\" annotation in metadata file conflicts with metadata convention. Convention expects \"group\" values.',
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        # invalid array-based metadata type: disease__time_since_onset
        self.assertIn(
            "disease__time_since_onset: 'three' in '36|three|1' does not match expected 'number' type.",
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        self.assertIn(
            "disease__time_since_onset: 'zero' in 'zero' does not match expected 'number' type.",
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        # invalid boolean value: disease__treated
        self.assertIn(
            "disease__treated: 'T' in 'T|F' does not match expected 'boolean' type.",
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        self.assertIn(
            "disease__treated: 'F' in 'F' does not match expected 'boolean' type.",
            metadata.issues["error"]["content"].keys(),
            "metadata validation should fail if metadata type does not match convention-designated type",
        )
        # non-uniform unit values: organism_age__unit
        self.assertIn(
            "disease__time_since_onset__unit: values for each unit metadata required to be uniform.",
            metadata.issues["error"]["content"].keys(),
            "Ontology_label values for units should be uniform",
        )
        # missing ontology ID or label for required metadata: organ
        self.assertIn(
            "organ: '' does not match '^[-A-Za-z0-9]+[_:][-A-Za-z0-9]+'",
            metadata.issues["error"]["convention"].keys(),
            "Missing required metadata should fail validation",
        )
        # missing ontology ID with existing label invalid even for non-required metadata: ethnicity
        self.assertIn(
            "'ethnicity' is a dependency of 'ethnicity__ontology_label'",
            metadata.issues["error"]["convention"].keys(),
            "missing ontology ID with existing label invalid even for non-required metadata",
        )
        # invalid header content: donor info (only alphanumeric or underscore allowed)
        self.assertIn(
            "donor info: only alphanumeric characters and underscore allowed in metadata name.",
            metadata.issues["error"]["format"].keys(),
            "metadata names must follow header rules (only alphanumeric or underscore allowed)",
        )
        self.teardown_metadata(metadata)

        # Arrays have NA values
        metadata = set_up_test("has_na_in_array.tsv")
        self.assertIn(
            "disease__time_since_onset: 'None' in 'None' does not match expected 'number' type.",
            metadata.issues["error"]["content"].keys(),
            "Non-numeric 'None' provided instead of numeric array should fail",
        )
        self.assertIn(
            "disease__treated: 'N/A' in 'True|N/A|False' does not match expected 'boolean' type.",
            metadata.issues["error"]["content"].keys(),
            "Non-boolean 'N/A' provided in boolean array should fail",
        )
        self.assertIn(
            "disease__treated: 'None' in 'FALSE|None' does not match expected 'boolean' type.",
            metadata.issues["error"]["content"].keys(),
            "Non-boolean 'None' provided in boolean array should fail",
        )
        self.teardown_metadata(metadata)

        # File has NA values in required fields
        metadata = set_up_test("has_na_in_required_fields.csv")
        # missing ontology ID or label for required metadata
        self.assertIn(
            "organ: required column \"organ\" missing data",
            metadata.issues["error"]["content"].keys(),
            "Missing required metadata should fail validation",
        )
        self.assertIn(
            "species: required column \"species__ontology_label\" missing data",
            metadata.issues["error"]["content"].keys(),
            "Missing required metadata should fail validation",
        )
        # NA value in required field (ontology)
        self.assertIn(
            "organ: '' does not match '^[-A-Za-z0-9]+[_:][-A-Za-z0-9]+'",
            metadata.issues["error"]["convention"].keys(),
            "Providing NaN for required metadata should fail validation",
        )
        # NA value in required field (enum)
        self.assertIn(
            "sex: 'Nan' is not one of ['male', 'female', 'mixed', 'unknown']",
            metadata.issues["error"]["convention"].keys(),
            "Providing NaN for required metadata should fail validation",
        )
        # NA value in non-required ontologyID field where ontology label provided
        self.assertIn(
            "'ethnicity' is a dependency of 'ethnicity__ontology_label'",
            metadata.issues["error"]["convention"].keys(),
            "Providing NaN for non-required ontologyID field where ontology label provided should fail validation",
        )
        self.teardown_metadata(metadata)

    def test_bigquery_json_content(self):
        """generated newline delimited JSON for BigQuery upload should match expected output
        """
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_array_v2.1.2.txt"
        )
        metadata, convention = self.setup_metadata(args)
        metadata.preprocess(is_metadata_convention=True)
        validate_input_metadata(metadata, convention, bq_json=True)

        generated_bq_json = str(metadata.study_file_id) + ".json"
        # This reference file needs updating with every new metadata convention version
        reference_bq_json = "../tests/data/bq_test.json"
        self.assertListEqual(
            list(io.open(generated_bq_json)), list(io.open(reference_bq_json))
        )

        self.teardown_metadata(metadata)

        # clean up downloaded generated BigQuery upload file
        try:
            os.remove("addedfeed000000000000000.json")
        except OSError:
            print("no file to remove")

    def test_invalid_mba_content(self):
        """Mouse Brain Atlas metadata should validate against MBA ontology file
        """
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/invalid_mba_v2.1.2.tsv"
        )
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
        #   missing organ_region when organ_region__ontology_label provided
        self.assertIn(
            "'organ_region' is a dependency of 'organ_region__ontology_label'",
            metadata.issues["error"]["convention"].keys(),
            "missing organ_region when organ_region__ontology_label provided should cause error",
        )
        #   Invalid identifier MBA_999999999
        self.assertIn(
            "organ_region: No match found in Allen Mouse Brain Atlas for provided ontology ID: MBA_999999999",
            metadata.issues["error"]["ontology"].keys(),
            "Invalid identifier MBA_999999999 should cause error",
        )
        #   mismatch of organ_region__ontology_label value with label value in MBA
        self.assertIn(
            'organ_region: input ontology_label \"Crus 1, urkinje layer\" does not match Allen Mouse Brain Atlas lookup \"Crus 1, Purkinje layer\" for ontology id \"MBA_000010676\".',
            metadata.issues["error"]["ontology"].keys(),
            "mismatch of organ_region__ontology_label value with label value in MBA should error",
        )
        #   mismatch of organ_region__ontology_label value with label from MBA_id lookup
        self.assertIn(
            'organ_region: input ontology_label \"Paraflocculus, granular layer\" does not match Allen Mouse Brain Atlas lookup \"Copula pyramidis, molecular layer\" for ontology id \"MBA_000010686\".',
            metadata.issues["error"]["ontology"].keys(),
            "mismatch of organ_region__ontology_label value with label from MBA_id lookup should error",
        )
        self.teardown_metadata(metadata)

    def test_is_empty_string(self):
        self.assertTrue(is_empty_string(""))
        self.assertTrue(is_empty_string("  "))
        self.assertFalse(is_empty_string("Hello"))
        self.assertFalse(is_empty_string(4))

    @patch("requests.get", side_effect=mocked_requests_get)
    def test_request_backoff_handling(self, mocked_requests_get):
        """errors in retrieving data from external resources should attempt MAX_HTTP_ATTEMPTS times and throw an exception
        """
        request_url = "https://www.ebi.ac.uk/ols/api/ontologies/"
        self.assertRaises(
            requests.exceptions.RequestException, request_json_with_backoff, request_url
        )
        self.assertEqual(mocked_requests_get.call_count, MAX_HTTP_ATTEMPTS)

    def test_is_label_or_synonym(self):
        label = "10x 3' v2"
        possible_matches = {
            "label": "10x 3' v2",
            "synonyms": ["10X 3' v2", "10x 3' v2 sequencing"],
        }
        self.assertTrue(is_label_or_synonym(possible_matches, label))
        label = "10X 3' v2"
        self.assertTrue(is_label_or_synonym(possible_matches, label))
        label = "10X 3' v2 sequencing"
        self.assertTrue(is_label_or_synonym(possible_matches, label))
        label = "10x 5' v3"
        self.assertFalse(is_label_or_synonym(possible_matches, label))

    def test_will_allow_synonym_matches(self):
        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/valid_no_array_synonyms_v2.0.0.txt"
        )
        metadata, convention = self.setup_metadata(args)
        self.assertTrue(
            metadata.validate_format(), "Valid metadata headers should not elicit error"
        )
        validate_input_metadata(metadata, convention)
        self.assertFalse(
            report_issues(metadata), "Valid ontology content should not elicit error"
        )
        self.teardown_metadata(metadata)

    def test_array_synonym_replacement(self):
        data = "../tests/data/annotation/metadata/convention/df.json"
        df = pd.read_json(data, lines=True)

        metadata_name = "ethnicity__ontology_label"
        matches_before_replace = [v == ["white"] for v in df[metadata_name]]

        replace_single_value_array(df, metadata_name, "white", "European")
        matches_after_replace = [v == ["white"] for v in df[metadata_name]]

        self.assertTrue(
            np.count_nonzero(matches_before_replace) == 1,
            "original df should have one instance of ['white']",
        )
        self.assertTrue(
            np.count_nonzero(matches_after_replace) == 0,
            "resulting df should have no instances of ['white']",
        )

        metadata_name = "disease__ontology_label"
        orig_values = list(df[metadata_name].transform(tuple).unique())
        replace = {"diabetes": "diabetes mellitus", "breast infection": "mastitis"}
        replace_synonym_in_multivalue_array(df, metadata_name, replace)
        replaced_values = list(df[metadata_name].transform(tuple).unique())

        expected_result = [
            ('diabetes mellitus', 'mastitis'),
            ('common cold',),
            ('diabetes mellitus', 'common cold'),
            ('diabetes mellitus', 'mastitis', 'common cold'),
            ('disease or disorder',),
        ]

        self.assertFalse(
            orig_values == replaced_values,
            "multi-value array names should be different after replacement",
        )

        self.assertEqual(
            replaced_values,
            expected_result,
            "multi-value array names should match expected result",
        )

    def test_validate_nonconventional_numeric_content(self):
        """Nonconventional numeric metadata values should all validate as numeric
        """

        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/invalid_nonconventional_numeric_v2.2.0.txt"
        )
        metadata, convention = self.setup_metadata(args)
        validate_input_metadata(metadata, convention)

        self.assertTrue(
            report_issues(metadata),
            "Numeric metadata with non-numeric value should fail validation",
        )

        self.assertEqual(
            list(metadata.issues["error"]["content"].keys())[0],
            'percent_mt: supplied value 07.juil is not numeric.',
            "expected error message not generated",
        )
        self.teardown_metadata(metadata)

    def test_excel_drag_check(self):
        """Evidence of an "Excel drag" event should generate errors
           Successful detection avoids EBI OLS queries when input data is faulty
        """

        args = (
            "--convention ../schema/alexandria_convention/alexandria_convention_schema.json "
            "../tests/data/annotation/metadata/convention/invalid_excel_drag.txt"
        )
        # Note: this test metadata file does not have array-based metadata
        # The file has optional ontology metadata lacking ontology labels.
        # Lack of labels triggers EBI OLS queries to "fill in" missing data.
        # This test will be slow due to those queries
        metadata, convention = self.setup_metadata(args)
        validate_input_metadata(metadata, convention)

        # The threshold for "excel drag" is currently >25 incrementing ontologyIDs.
        # The following two cases in the test file are not (yet) expected to fail:
        # - ethnicity has 25 consecutive incrementing ontologyIDs (all labeled "European")
        # - geographical_region has 25 consecutive incrementing ontologyIDs (unlabeled)
        # Note: all 25 of the geographical_region incrementing IDs are valid ontologyIDs

        ethnicity_error = False
        geographical_region_error = False
        for error in metadata.issues["error"]["ontology"].keys():
            if "ethnicity" in error:
                ethnicity_error = True
            if "geographical_region_error" in error:
                geographical_region_error = True
        self.assertFalse(
            ethnicity_error,
            "25 consecutive incrementing ontologyIDs should not trigger error",
        )
        self.assertFalse(
            geographical_region_error,
            "25 consecutive incrementing ontologyIDs should not trigger error",
        )

        # Metadata 'race' has multiple ontologyIDs paired with the same label,
        # note: all incrementing ontologyIDs for race are valid IDS AND the
        # input data has no ontology_label for race
        # For most SCP-required ontologies, a casual search
        # did not find runs of >25 actual, valid ontologyIDs
        # It seems reasonable to mark 25 consecutive, adjacent ontologyIDs
        # in cell-based data as error (ie. highly unlikely by chance)
        self.assertIn(
            "race: Long stretch of contiguously incrementing ontology ID values suggest cut and paste issue - "
            "exiting validation, ontology content not validated against ontology server.\n"
            "Please confirm ontology IDs are correct and resubmit.\n",
            metadata.issues["error"]["ontology"].keys(),
            "Run of >25 incrementing ontology labels should fail validation",
        )

        # Metadata 'disease' also has multiple ontologyIDs paired with the same label
        # It advises to only check mismatches with labels that are truly multiply assigned:
        # ['disease or disorder', 'absent']
        self.assertIn(
            "disease: Long stretch of contiguously incrementing ontology ID values suggest cut and paste issue - "
            "exiting validation, ontology content not validated against ontology server.\n"
            "Please confirm ontology IDs are correct and resubmit.\n"
            "Check for mismatches between ontology ID and provided ontology label(s) "
            "['absent', 'disease or disorder'].\n",
            metadata.issues["error"]["ontology"].keys(),
            "ontology label multiply paired with IDs should error",
        )

        # Metadata 'species' has multiple ontologyIDs paired with the same label,
        # this triggers the additional statement to check for mismatches
        species_error = False
        for error in metadata.issues["error"]["ontology"].keys():
            if "species" in error:
                species_error = True
        self.assertTrue(
            species_error,
            "ontology label multiply paired with IDs should error (depends on EBI OLS availability)",
        )

        self.teardown_metadata(metadata)


if __name__ == "__main__":
    unittest.main()
