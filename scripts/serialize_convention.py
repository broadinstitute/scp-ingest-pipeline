"""Produce JSON Schema from metadata convention tsv file

DESCRIPTION
This CLI takes a tsv metadata convention and creates a JSON Schema representation.
The JSON schema represents the rules that should be enforced on metadata files
for studies participating under the convention.
This CLI also generates a JSON file representing the BigQuery schema and incorporates
info specifying the non-convention, SCP-added BigQuery terms (scp_bq_inputs.json)

EXAMPLE
# Generate JSON file for Alexandria from the Alexandria metadata convention TSV
# Expects TSV and scp_bq_inputs.json files in location under
# scp-ingest-pipeline/schema/<project>_convention/snapshot/<version>
$ python serialize_convention.py alexandria 2.0.0

"""

import argparse
import csv
import os
import json
import re
import requests

from convention_to_bq_schema import (
    load_scp_inputs,
    build_bq_schema,
    add_scp_fields_to_schema,
    write_bq_schema,
)


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention TSV file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # add arguments
    parser.add_argument(
        'project', help='One-word project name that metadata convention belongs to'
    )
    parser.add_argument('version', help='semantic version number for convention')
    return parser


def add_dependency(key, value, dict):
    """
    Add dependency to appropriate dictionary

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
    # handle controlled-list attributes as enum
    if row['class'] == 'enum':
        dict['type'] = row['type']
        dict[row['class']] = row['controlled_list_entries']
    elif row['class'] == 'ontology':
        dict['type'] = row['type']
        ontology_format = r'^[-A-Za-z0-9]+[_:][-A-Za-z0-9]+'
        dict['pattern'] = ontology_format
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

    # include syntax constraint for ontologies
    # ontologies whose short names include underscore or colon will fail to validate
    elif row['class'] == 'ontology':
        ontology_format = r'^[A-Za-z]+[_:][0-9]'
        dict['pattern'] = ontology_format
        dict['type'] = row['type']

    else:
        dict['type'] = row['type']


def build_schema_info(project, version):
    """
    generate dictionary of schema info for the project
    """
    info = {}
    info['$schema'] = 'https://json-schema.org/draft-07/schema#'
    # $id below is a placeholder, not functional yet
    info['$id'] = (
        f'https://singlecell.broadinstitute.org/single_cell/api/v1/metadata-schema/'
        f'{project}_convention/{version}/{project}_convention_schema.json'
    )
    info['title'] = project + ' metadata convention'
    info['description'] = 'metadata convention for the ' '%s project' % (project)
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
        # Values for enums are in an array and json.dump adds quotes around
        # the array representation and escapes the internal quotes
        # which both cause invalid JSON schema
        jsonstring = re.sub(r'"\[', r'[', jsonstring)
        jsonstring = re.sub(r'"\]"', r'"]', jsonstring)
        jsonstring = re.sub(r'\\"', r'"', jsonstring)
    return jsonstring


def write_json_schema(filename, object):
    """
    write JSON Schema file
    """
    with open(filename, 'w') as jsonfile:
        jsonfile.write(object)


def set_file_names(project, version):
    """
    Infer input tsv and scp_bq_inputs.json file locations from project, version info
    If input tsv doesn't exist, exit
    If scp_bq_inputs.json doesn't exist, exit
    If metadata convention file already exists, exit
    """
    inputname = f'../schema/{project}_convention/snapshot/{version}/{project}_convention_schema.tsv'
    convention_filename = f'../schema/{project}_convention/snapshot/{version}/{project}_convention_schema.json'
    bq_schema_filename = f'../schema/{project}_convention/snapshot/{version}/{project}_convention_schema.bq_schema.json'
    scp_bq_inputs = (
        f'../schema/{project}_convention/snapshot/{version}/scp_bq_inputs.json'
    )
    if not os.path.exists(inputname):
        print(f'{inputname} does not exist, please check your TSV file and try again')
        exit(1)
    if not os.path.exists(scp_bq_inputs):
        print(
            f'{scp_bq_inputs} does not exist, please check your JSON file and try again'
        )
        exit(1)
    if os.path.exists(convention_filename):
        print(f'{convention_filename} already exists, please delete file and try again')
        exit(1)
    if os.path.exists(bq_schema_filename):
        print(f'{bq_schema_filename} already exists, please delete file and try again')
        exit(1)
    filenames = {}
    filenames['input'] = inputname
    filenames['convention_filename'] = convention_filename
    filenames['scp_bq_inputs'] = scp_bq_inputs
    return filenames


def retrieve_ontology(ontology_url):
    """Retrieve an ontology listing from EBI OLS
    :param ontology_term: identifier of a term in an ontology in OLS (e.g. CL_0002419)
    :return: JSON payload of ontology, or None
    """
    response = requests.get(ontology_url)
    if response.status_code == 200:
        return response.json()
    else:
        return None


def serialize_convention(convention, input_tsv):
    """
    Build convention as a Python dictionary
    """

    properties = {}
    required = []
    dependencies = {}

    with open(input_tsv) as tsvfile:
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
                add_dependency(row['attribute'], row['dependency'], dependencies)
            if row['dependent']:
                # dependent is bi-directional (aka "required-if")
                add_dependency(row['attribute'], row['dependent'], dependencies)
                add_dependency(row['dependent'], row['attribute'], dependencies)

            # build dictionary for each attribute
            if row['attribute_description']:
                entry['description'] = row['attribute_description']

            # handle properties unique to the ontology class of attributes
            # must allow for multiple ontologies, comma-delimited
            if row['class'] == 'ontology':
                ontology_urls = row['ontology'].split(',')
                if len(ontology_urls) > 1:
                    print(f"Checking {len(ontology_urls)} urls for {row['attribute']}")
                for url in ontology_urls:
                    ontology = retrieve_ontology(url)
                    if ontology:
                        print("Valid ontology URL", url)
                    else:
                        print("Invalid ontology URL", url)
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


def write_convention_schema(dict, filepath):
    dump_json(dict, filepath)
    write_json_schema(filepath, clean_json(filepath))


def generate_bq_schema(filenames):
    scp_bq_input = load_scp_inputs(filenames['scp_bq_inputs'])
    schema = build_bq_schema(filenames['input'])
    schema = add_scp_fields_to_schema(schema, scp_bq_input)
    write_bq_schema(schema, filenames['convention_filename'])


if __name__ == '__main__':
    args = create_parser().parse_args()
    project = args.project
    version = args.version
    filenames = set_file_names(project, version)
    schema_info = build_schema_info(project, version)
    convention = serialize_convention(schema_info, filenames['input'])
    write_convention_schema(convention, filenames['convention_filename'])
    generate_bq_schema(filenames)
    print(f"Metadata convention written to {filenames['convention_filename']}")
