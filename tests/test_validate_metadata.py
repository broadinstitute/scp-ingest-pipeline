import sys
import unittest

sys.path.append('../ingest')
sys.path.append('../ingest/validation')

from validate_metadata import *


class TestValidateMetadata(unittest.TestCase):
    def setup_metadata(self, args):
        args_list = args.split(' ')
        args = create_parser().parse_args(args_list)
        with open(args.convention, 'r') as f:
            convention = json.load(f)
        filetsv = args.input_metadata
        metadata = CellMetadata(filetsv)
        metadata.validate_format()
        return (metadata, convention)

    def test_header_format(self):
        """Header rows of metadata file should conform to standard
        """

        args = (
            '../tests/data/AMC_v0.8.json '
            '../tests/data/error_headers.tsv'
        )
        metadata = self.setup_metadata(args)[0]
        self.assertFalse(metadata.validate_header_keyword())
        self.assertIn(
            'Error: Metadata file header row malformed, missing NAME',
            metadata.errors['format'],
            'Validation should detect malformed NAME keyword'
        )

        self.assertFalse(metadata.validate_type_keyword())
        self.assertIn(
            'Error:  Metadata file TYPE row malformed, missing TYPE',
            metadata.errors['format'],
            'Validation should detect malformed TYPE keyword'
        )

        self.assertFalse(
            metadata.validate_type_annotations(),
            'Validation should restrict type annotations to allowed values'
        )

    def test_convention_content(self):
        """Metadata convention should be valid jsonschema
            """

        args = (
            '../tests/data/AMC_invalid.json '
            '../tests/data/metadata_valid.tsv'
        )
        convention = self.setup_metadata(args)[1]
        self.assertIsNone(
            validate_schema(convention),
            'Invalid metadata schema should be detected'
        )

    def test_nonontology_content(self):
        """Non-ontology metadata should conform to convention requirements
            """
        args = (
            '../tests/data/AMC_v0.8.json '
            '../tests/data/metadata_valid.tsv'
        )
        metadata, convention = self.setup_metadata(args)
        metadata_valid = metadata.validate_format()
        self.assertTrue(
            metadata.validate_format(),
            'Valid metadata headers should not elicit error'
        )
        if metadata_valid:
            process_metadata_content(metadata, convention)
        self.assertFalse(
            report_errors(metadata),
            'Valid metadata content should not elicit error'
        )


if __name__ == '__main__':
    unittest.main()
