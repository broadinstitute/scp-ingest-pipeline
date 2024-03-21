"""Produce BigQuery schema JSON file from metadata convention tsv file

DESCRIPTION
This CLI takes a tsv metadata convention and creates a BigQuery schema JSON file.
JSON file is written to the directory where the input file lives.

SYNTAX
$ python convention_to_bq_schema.py scp_bq_inputs.json metadata_convention.tsv

EXAMPLE
$ python convention_to_bq_schema.py ../schema/alexandria_convention/snapshot/2.0.0/scp_bq_inputs.json ../schema/alexandria_convention/snapshot/2.0.0/alexandria_convention_schema.tsv

"""

import argparse
import csv
import os
import json


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
    parser.add_argument(
        'scp_input', help='SCP requirements for BigQuery schema, JSON file'
    )
    parser.add_argument('input_convention', help='Metadata convention tsv file')
    return parser


def process_row_type(type_info):
    type_map = {
        'integer': 'integer',
        'boolean': 'bool',
        'string': 'string',
        'number': 'float',
    }
    return type_map.get(type_info).upper()


def build_bq_schema(input_convention):
    """
    Build schema as a Python dictionary
    """
    #
    with open(input_convention) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        schema = []
        for row in reader:
            entry = {}
            entry['name'] = row['attribute']
            entry['type'] = process_row_type(row['type'])
            # handle arrays of values
            if row['array']:
                entry['type'] = process_row_type(row['type'])
                entry['mode'] = 'REPEATED'
            else:
                entry['mode'] = 'NULLABLE'
            schema.append(entry)
    return schema


def add_scp_fields_to_schema(schema, scp_inputs):
    for entry in scp_inputs:
        schema_entry = {
            'name': entry,
            'type': scp_inputs[entry][0],
            'mode': scp_inputs[entry][1],
        }
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


def write_bq_schema(data, inputname, filepath=''):
    """
    write BigQuery schema as json file
    """
    filename = generate_output_name(inputname, filepath)
    with open(filename, 'w') as jsonfile:
        json.dump(data, jsonfile, sort_keys=True, indent=4)


def load_scp_inputs(inputfile):
    with open(inputfile, 'r') as f:
        scp_inputs = json.load(f)
    return scp_inputs


if __name__ == '__main__':
    args = create_parser().parse_args()
    scp_input = args.scp_input
    scp_bq_input = load_scp_inputs(scp_input)
    input_convention = args.input_convention
    output_path = args.output_path
    schema = build_bq_schema(input_convention)
    schema = add_scp_fields_to_schema(schema, scp_bq_input)
    write_bq_schema(schema, input_convention, output_path)
