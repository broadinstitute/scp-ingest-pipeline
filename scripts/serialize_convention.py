#! /usr/bin/python

import argparse
import csv
import os
import json
import re

"""Produce JSON Schema from metadata convention tsv file

DESCRIPTION
This CLI takes a tsv metadata convention and creates a JSON Schema representation.
The JSON schema represents the rules that should be enforced on metadata files 
for studies participating under the convention.

EXAMPLE
# Generate json file for Alexandria from the Alexandria metadata convention tsv
$ serialize_convention.py project metadata_convention.tsv 

# if using the test data in this repo and wanted output in current working dir
$ serialize_convention.py -p . Alexandria ../tests/data/AMC_v0.8.tsv 

"""


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention tsv file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # add arguments
    parser.add_argument('--label', '-l', type=str,
                        help='Label to insert into the file name ' +
                        'for the Metadata convention json file [optional]',
                        default=None),
    parser.add_argument('--output-path', '-p', type=str,
                        help='Path for output file [optional]', default=None)
    parser.add_argument('project', type=str,
                        help='One-word project name that the metadata ' +
                        'convention belongs to [Required]')
    parser.add_argument('input_convention',
                        help='Metadata convention tsv file [Required]')
    return parser


def add_dependency(key, value, dict):
    """
    Add dependency to appopriate dictionary

    ToDo: check if defaultdict would eliminate this function
    """
    if key in dict:
        if value not in dict[key]:
            dict[key].append(value)
    else:
        dict[key] = [value]


def build_array_object(row):
    """
    Build "items" dictionary object according to Class type
    """
    dict = {}
    # handle time attributes as string with time format
    if row['class'] == 'time':
        dict['type'] = row['type']
        dict['format'] = row['class']
        return dict
    # handle controlled-list attributes as enum
    elif row['class'] == 'enum':
        dict['type'] = row['type']
        dict[row['class']] = row['controlled_list_entries']
        return dict
    else:
        dict['type'] = row['type']
        return dict


def build_single_object(row, dict):
    """
    Add appropriate class properties for non-array attributes
    """
    # handle time attributes as string with time format
    if row['class'] == 'time':
        dict['type'] = row['type']
        dict['format'] = row['class']

    # handle controlled-list attributes as enum
    elif row['class'] == 'enum':
        dict['type'] = row['type']
        dict[row['class']] = row['controlled_list_entries']

    else:
        dict['type'] = row['type']


def build_schema_info(project):
    """
    generate dictionary of schema info for the project
    """
    info = {}
    info['$schema'] = 'https://json-schema.org/draft-07/schema#'
    # $id below is a placeholder, not functional yet
    info['$id'] = ('https://singlecell.broadinstitute.org/api/v1/metadata-schemas/'
                   '%s.schema.json' % (project))
    info['title'] = project + ' Metadata Convention'
    info['description'] = ('Metadata convention for the '
                           '%s project' % (project))
    return info


def dump_json(dict, filename):
    """
    write metadata convention json file
    """
    with open(filename, 'w') as jsonfile:
        json.dump(dict, jsonfile, sort_keys=True, indent=4)


def clean_json(filename):
    """
    remove escape characters to produce proper JSON Schema format
    """
    with open(filename, 'r') as jsonfile:
        jsonstring = jsonfile.read()
        jsonstring = re.sub(r'"\[', r'[', jsonstring)
        jsonstring = re.sub(r'\]"', r']', jsonstring)
        jsonstring = re.sub(r'\\"', r'"', jsonstring)
    return jsonstring


def write_json_schema(filename, object):
    """
    write JSON Schema file
    """
    with open(filename, 'w') as jsonfile:
        jsonfile.write(object)


def generate_output_name(inputname, label, path = ''):
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


def serialize_convention(convention, input_convention):
    """
    Build convention as a Python dictionary
    """

    properties = {}
    required = []
    dependencies = {}

    with open(input_convention) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')

        for row in reader:
            entry = {}

            # build list of required attributes for metadata convention schema
            if row['required']:
                required.append(row['attribute'])

            # build dictionary of dependencies for metadata convention schema
            if row['dependency']:
                # dependencies (aka "if" relationships) are uni-directional
                # if 'attribute', 'dependency' must also exist
                add_dependency(
                    row['attribute'], row['dependency'], dependencies)
            if row['dependent']:
                # dependent is bi-directional (aka "required-if")
                add_dependency(row['attribute'], row['dependent'], dependencies)
                add_dependency(row['dependent'], row['attribute'], dependencies)

            # build dictionary for each attribute
            if row['default']:
                entry['default'] = row['default']
            if row['attribute_description']:
                entry['description'] = row['attribute_description']

            # handle properties unique to the ontology class of attributes
            if row['class'] == 'ontology':
                entry[row['class']] = row['ontology']

            # handle arrays of values
            if row['array']:
                entry['type'] = 'array'
                entry['items'] = build_array_object(row)
            else:
                build_single_object(row, entry)

            if row['dependency_condition']:
                entry['dependency_condition'] = row['dependency_condition']

            # build dictionary of properties for the metadata convention schema
            properties[row['attribute']] = entry

        # build metadata convention schema
        convention['properties'] = properties
        convention['required'] = required
        convention['dependencies'] = dependencies
    return convention


def write_schema(dict, inputname, label, filepath = ''):
    filename = generate_output_name(inputname, label, filepath)
    dump_json(dict, filename)
    write_json_schema(filename, clean_json(filename))


if __name__ == '__main__':
    args = create_parser().parse_args()
    input_convention = args.input_convention
    label = args.label
    output_path = args.output_path
    project = args.project
    schema_info = build_schema_info(project)
    convention = serialize_convention(schema_info, input_convention)
    write_schema(convention, input_convention, label, output_path)
