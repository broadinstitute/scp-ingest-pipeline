#! /usr/bin/python
"""Validate input metadata tsv file against metadata convention

DESCRIPTION
This CLI takes a tsv metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using json file for Alexandria metadata convention tsv, validate input tsv
$ validate_metadata.py AMC_v0.8.json metadata_test.tsv

"""

import csv
import argparse
import re
import json
import jsonschema
from collections import defaultdict


class Cell_Metadata:
    def __init__(self, file_path):
        self.file = open(file_path, 'r')
        self.headers = self.file.readline().rstrip('\n').split("\t")
        self.metadata_types = self.file.readline().rstrip('\n').split("\t")
        self.annotation_type = ['group', 'numeric']
        self.errors = defaultdict(list)

    def validate_header_keyword(self):
        if self.headers[0].casefold() == 'NAME'.casefold():
            try:
                self.headers.pop(self.headers.index('NAME'))
            except ValueError:
                print(
                    'Warning: metadata file keyword "NAME" provided as {x}'.
                    format(x=self.headers[0])
                )
                self.headers.pop(0)
        else:
            self.errors['format'].append(
                ('Error: Metadata file header row malformed, missing NAME')
            )
        return

    def validate_type_keyword(self):
        if self.metadata_types[0].casefold() == 'TYPE'.casefold():
            try:
                self.metadata_types.pop(self.metadata_types.index('TYPE'))
            except ValueError:
                print(
                    'Warning: metadata file keyword "TYPE" provided as {x}'.
                    format(x=self.metadata_types[0])
                )
                self.metadata_types.pop(0)
        else:
            self.errors['format'].append(
                ('Error:  Metadata file TYPE row malformed, missing TYPE')
            )
        return

    def validate_type_annotations(self):
        annot_err = False
        annots = []
        for t in self.metadata_types:
            if t not in self.annotation_type:
                annots.append(t)
                annot_err = True
        if annot_err:
            self.errors['format'].append(
                (
                    'Error: TYPE declarations should be "group" or "numeric"; '
                    'Please correct: {annots}'.format(
                        annots=', '.join(map(str, annots))
                    )
                )
            )
        return

    def validate_against_header_count(self, list):
        if not len(self.headers) == len(list):
            self.errors['format'].append(
                str(
                    'Error: {x} TYPE declarations for {y} column headers'.
                    format(x=len(self.headers), y=len(list))
                )
            )
        return

    def validate_format(self):
        self.validate_header_keyword()
        self.validate_type_keyword()
        self.validate_type_annotations()
        self.validate_against_header_count(self.metadata_types)
        if self.errors['format']:
            valid = False
        else:
            valid = True
        return valid


def create_parser():
    """
    Parse command line values for validate_metadata

    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # parser.add_argument(
    #     '--output',
    #     '-o',
    #     type=str,
    #     help='Output file name [optional]',
    #     default=None
    # )
    parser.add_argument(
        'convention', type=str, help='Metadata convention json file [Required]'
    )
    parser.add_argument('input_metadata', help='Metadata tsv file [Required]')
    return parser


def load_schema(schemafile):
    """
    Read Convention
    :param schemafile: metadata convention file
    :return: dict representing metadata convention
    """
    with open(schemafile, "r") as read_file:
        schema = json.load(read_file)


# ToDo - action if schema is invalid?
    jsonschema.Draft7Validator.check_schema(schema)
    valid_schema = jsonschema.Draft7Validator(schema)

    return valid_schema
"""
Read tsv input row by row
"""
"""
ontology validation
"""
"""
generate error report
"""
"""
WAIT: handle array data types

"""
"""
DEFER (loom): Check loom format is valid
  what are the criteria?
"""
"""
DEFER (loom): Read loom metadata, row by row?
"""
"""
DEFER: Things to check before pass intended data types to FireStore
ensure numeric (TYPE == numeric in tsv; type number or integer in loom)
    stored as numeric
ensure group (even if it is a number) stored as string
NaN stored as null?
error on empty cells
"""
"""

"""

if __name__ == '__main__':
    args = create_parser().parse_args()
    arguments = vars(args)
    schema = load_schema(args.convention)
    filetsv = args.input_metadata
    metadata = Cell_Metadata(filetsv)
    print('Validating', filetsv)
    metadata.validate_format()
    print("error:", metadata.errors.items())

    # compiled regex to identify arrays in metadata file
    array_format = re.compile(r'\[.*\]')

    with open(filetsv) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        # skip TYPE row
        type = next(reader)
        for row in reader:
            #    print(row['NAME'], row['disease'])
            # DictReader values are strings
            #   reformat all intended arrays from string
            # TODO replace empty values with NONE?
            #   or other handling for empties
            for key, value in row.items():
                if not value:
                    row[key] = None
                if array_format.match(value):
                    row[key] = json.loads(value)
                if type[key] == "numeric":
                    try:
                        row[key] = float(value)
                    # terrible antipattern below - must fix
                    except Exception:
                        pass

            for error in schema.iter_errors(row):
                print(row['NAME'], "error:", error.message)
#        print(row)
