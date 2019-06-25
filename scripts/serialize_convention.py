#! /usr/bin/python

import argparse
import csv
import os
import json
import re

"""Produce JSON Schema from metadata convention tsv file"""

__author__ = "Jean Chang"
__copyright__ = "Copyright 2019"
__license__ = "MIT"
__email__ = "jlchang@broadinstitute.org"
__status__ = "Development"


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention tsv file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description='Produce JSON Schema from metadata convention tsv file.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # add arguments
    parser.add_argument('input_convention',
                        help='Metadata convention tsv file')
    parser.add_argument('--collective', '-c', type=str,
                        help='name of the project that the metadata ' +
                        'convention belongs to [Required]',
                        required=True, default=None)
    parser.add_argument('--label', '-l', type=str,
                        help='Label to insert into the file name ' +
                        'for the Metadata convention json file [optional]',
                        default=None),
    parser.add_argument('--output_file', '-o', type=str,
                        help='Output file name [optional]', default=None)
    return parser


def addDependency(key, value, dict):
    """
    Add dependency to appopriate dictionary

    ToDo: check if defaultdict would eliminate this function
    """
    if key in dict:
        if value not in dict[key]:
            dict[key].append(value)
    else:
        dict[key] = [value]


def buildArrayObject(row):
    """
    Build "items" dictionary object according to Class type
    """
    dict = {}
    # handle time attributes as string with time format
    if row['class'] == 'time':
        dict['type'] = row['type']
        dict['format'] = row['class']
        print('dict time', dict, 'for', row['attribute'])
        return dict
    # handle controlled-list attributes as enum
    elif row['class'] == 'enum':
        dict['type'] = row['type']
        dict[row['class']] = row['controlled_list_entries']
        print('dict enum', dict, 'for', row['attribute'])
        return dict
    else:
        dict['type'] = row['type']
        return dict


def buildSingleObject(row, dict):
    """
    Add appropriate class properties for non-array attributes
    """
    # handle time attributes as string with time format
    if row['class'] == 'time':
        dict['type'] = row['type']
        dict['format'] = row['class']
        print('single time', dict, 'for', row['attribute'])

    # handle controlled-list attributes as enum
    elif row['class'] == 'enum':
        dict['type'] = row['type']
        dict[row['class']] = row['controlled_list_entries']
        print('single enum', dict, 'for', row['attribute'])

    else:
        dict['type'] = row['type']


def buildSchemaInfo(collective):
    """
    generate dictionary of schema info for the collective
    """
    info = {}
    info['$schema'] = 'http://json-schema.org/draft-07/schema#'
    info['$id'] = ('http://singlecell.broadinstitute.org/schemas/'
                   '%s.schema.json' % (collective))
    info['title'] = collective + ' Metadata Convention'
    info['description'] = ('Metadata convention for the '
                           '%s project' % (collective))
    return info


def dumpJSON(dict, filename):
    """
    write metadata convention json file
    """
    with open(filename, 'w') as jsonfile:
        json.dump(dict, jsonfile, sort_keys=True, indent=4)
    jsonfile.close()
    print("end dumpJSON")


def cleanJSON(filename):
    """
    remove escape characters to produce proper JSON Schema format
    """
    with open(filename, 'r') as jsonfile:
        jsonstring = jsonfile.read()
        jsonstring = re.sub(r'"\[', r'[', jsonstring)
        jsonstring = re.sub(r'\]"', r']', jsonstring)
        jsonstring = re.sub(r'\\"', r'"', jsonstring)
    jsonfile.close()
    return jsonstring


def writeJSONSchema(filename, object):
    """
    write JSON Schema file
    """
    with open(filename, 'w') as jsonfile:
        jsonfile.write(object)
    print("end writeJSONSchema")


def generateOutputName(inputname, label):
    """
    Build output filename from inputname
    """
    head, tail = os.path.split(inputname)
    name, suffix = os.path.splitext(tail)
    if label:
        labeledName = '.'.join([name, label, 'json'])
    else:
        labeledName = '.'.join([name, 'json'])
    if head:
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
                addDependency(
                    row['attribute'], row['dependency'], dependencies)
            if row['dependent']:
                # dependent is bi-directional (aka "required-if")
                addDependency(row['attribute'], row['dependent'], dependencies)
                addDependency(row['dependent'], row['attribute'], dependencies)

            # build dictionary for each attribute
            if row['default']:
                entry['default'] = row['default']
            if row['attribute_description']:
                entry['description'] = row['attribute_description']

            # handle properties unique to the ontology class of attributes
            if row['class'] == 'ontology':
                entry[row['class']] = row['ontology']
                entry['ontology_root'] = row['ontology_root']

            # handle arrays of values
            if row['array']:
                entry['type'] = 'array'
                entry['items'] = buildArrayObject(row)
            else:
                buildSingleObject(row, entry)

            if row['dependency_condition']:
                entry['dependency_condition'] = row['dependency_condition']

            # build dictionary of properties for the metadata convention schema
            properties[row['attribute']] = entry

        # build metadata convention schema
        convention['properties'] = properties
        convention['required'] = required
        convention['dependencies'] = dependencies

    tsvfile.close()
    return convention


def writeSchema(dict, inputname, label, filename):
    if filename:
        filename = generateOutputName(filename, label)
    else:
        filename = generateOutputName(inputname, label)
    dumpJSON(dict, filename)
    writeJSONSchema(filename, cleanJSON(filename))


if __name__ == '__main__':
    args = create_parser().parse_args()
    input_convention = args.input_convention
    label = args.label
    output_file = args.output_file
    collective = args.collective
    schemaInfo = buildSchemaInfo(collective)
    convention = serialize_convention(schemaInfo, input_convention)
    writeSchema(convention, input_convention, label, output_file)
