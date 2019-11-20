"""Produce BigQuery schema JSON file from metadata convention tsv file

DESCRIPTION
This CLI takes a tsv metadata convention and creates a BigQuery schema JSON file.
JSON file is written to the directory where the input file lives.

SYNTAX
$ python convention_to_bq_schema.py metadata_convention.tsv

EXAMPLE
$ python convention_to_bq_schema.py ../../tests/data/AMC_v1.1.3.tsv

"""

import argparse
import csv
import os
import json

REQUIRED_FIELDS = ['CellID', 'study_accession', 'file_id']


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention tsv file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--output-path', '-p', help='Path for output file')
    parser.add_argument('input_convention', help='Metadata convention tsv file')
    return parser


def process_row_type(type_info):
    type_map = {
        'integer': 'integer',
        'boolean': 'boolean',
        'string': 'string',
        'number': 'float',
    }
    return type_map.get(type_info).upper()


def build_schema(input_convention):
    """
    Build schema as a Python dictionary
    """

    with open(input_convention) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        schema = []
        for row in reader:
            entry = {}
            entry['name'] = row['attribute']
            entry['type'] = process_row_type(row['type'])
            if row['attribute'] in REQUIRED_FIELDS:
                entry['mode'] = 'REQUIRED'
            # handle arrays of values
            if row['array']:
                entry['type'] = process_row_type(row['type'])
                entry['mode'] = 'REPEATED'
            else:
                entry['mode'] = 'NULLABLE'
            schema.append(entry)
    return schema


def add_scp_fields_to_schema(schema):
    entries = ['study_accession', 'file_id']
    for entry in entries:
        schema_entry = {'name': entry, 'type': 'string', 'mode': 'REQUIRED'}
        schema.append(schema_entry)
    return schema


def generate_output_name(inputname, path='', label='bq_schema'):
    """
    Build output filename from inputname
    """
    head, tail = os.path.split(inputname)
    name, suffix = os.path.splitext(tail)
    if label:
        labeledName = '.'.join([name, label, 'json'])
    else:
        labeledName = '.'.join([name, 'json'])
    if path:
        outputname = '/'.join([path, labeledName])
    elif head:
        outputname = '/'.join([head, labeledName])
    else:
        outputname = labeledName
    return outputname


def write_schema(data, inputname, filepath=''):
    """
    write BigQuery schema as json file
    """
    filename = generate_output_name(inputname, filepath)
    with open(filename, 'w') as jsonfile:
        json.dump(data, jsonfile, sort_keys=True, indent=4)


if __name__ == '__main__':
    args = create_parser().parse_args()
    input_convention = args.input_convention
    output_path = args.output_path
    schema = build_schema(input_convention)
    schema = add_scp_fields_to_schema(schema)
    write_schema(schema, input_convention, output_path)
