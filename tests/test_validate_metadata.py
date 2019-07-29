import sys
import unittest

sys.path.append('../ingest')
sys.path.append('../ingest/validation')
from validate_metadata import *


class TestValidateMetadata(unittest.TestCase):
    def setup_validation(self, args):
        args_list = args.split(' ')
        args = create_parser().parse_args(args_list)
        filetsv = args.input_metadata
        metadata = Cell_Metadata(filetsv)
        metadata.validate_format()
        metadata.file.close()
        return metadata

    def test_format_name(self):
        """Header row of metadata file should start with NAME keyword
        """

        args = ('../tests/data/AMC_v0.8.json ' '../tests/data/error_test10.tsv')
        metadata = self.setup_validation(args)
        self.assertFalse(metadata.validate_header_keyword())

        self.assertIn(
            'Error: Metadata file header row malformed, missing NAME',
            metadata.errors['format'],
            "Validation should detect malformed NAME keyword"
        )

    def test_format_type(self):
        """Second row of metadata file should start with TYPE keyword
        """

        args = ('../tests/data/AMC_v0.8.json ' '../tests/data/error_test11.tsv')
        metadata = self.setup_validation(args)

        self.assertFalse(metadata.validate_type_keyword())
        self.assertIn(
            'Error:  Metadata file TYPE row malformed, missing TYPE',
            metadata.errors['format'],
            "Validation should detect malformed TYPE keyword"
        )

    def test_format_type_annotations(self):
        """TYPE values restricted to "group" or "numeric"
        """

        args = ('../tests/data/AMC_v0.8.json ' '../tests/data/error_test12.tsv')
        metadata = self.setup_validation(args)

        self.assertFalse(metadata.validate_type_annotations())


if __name__ == '__main__':
    unittest.main()
