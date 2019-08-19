#! /usr/bin/python
"""Validate input metadata TSV file against metadata convention.

DESCRIPTION
This CLI takes a TSV metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using JSON file for Alexandria metadata convention TSV, validate input TSV
$ validate_metadata.py AMC_v0.8.json metadata_test.tsv

"""

import argparse
import json
from collections import defaultdict
import logging
import sys

import jsonschema

sys.path.append('..')
from cell_metadata import *

# ToDo set up parameters to adjust log levels
#  logging.basicConfig(level=logging.INFO, filename='app.log', filemode='w',
#   format='%(name)s - %(levelname)s - %(message)s')
logging.basicConfig(
    level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_parser():
    """Parse command line values for validate_metadata
    """
    logger.debug('Begin: create_parser')
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
    parser.add_argument('convention', help='Metadata convention JSON file ')
    parser.add_argument('input_metadata', help='Metadata TSV file')
    return parser


def validate_schema(json):
    """Check validity of metadata convention as JSON schema.

    :param schemafile: metadata convention JSON file
    :return: if valid, jsonschema validator object using input convention
                or returns None
    """

    try:
        jsonschema.Draft7Validator.check_schema(json)
        valid_schema = jsonschema.Draft7Validator(json)
        return valid_schema
    except jsonschema.SchemaError as e:
        ### save this output as part of error output
        print('Input JSON is invalid as jsonschema;', e)
        return None


def extract_numeric_headers(metadata):
    """Find metadata headers of type numeric.

    ASSUMES metadata.validate_format() == True
    If the file headers are improperly formatted, extract_numeric_headers
    may not be able to correctly identify the headers intended to be numeric

    :param metadata: cell metadata object
    list of numeric-type metadata headers at metadata.type['numeric_headers']
    """
    logger.debug('Begin: extract_numeric_headers')
    for index, mtype in enumerate(metadata.metadata_types):
        if mtype == 'numeric':
            metadata.type['numeric_headers'].append(metadata.headers[index])
    return


def extract_convention_types(convention, metadata):
    """Populates metadata.type with property type from metadata convention

    :param convention: dict representation of metadata convention
    :param metadata: cell metadata object
    """
    logger.debug('Begin: extract_numeric_from_convention')
    # setup of metadata.type['convention']['ontology'] has onotology labels
    metadata.type['convention'] = defaultdict(list)
    for k in convention['properties'].keys():
        if convention['properties'][k]['type'] == 'number':
            metadata.type['convention']['number'].append(k)
        elif convention['properties'][k]['type'] == 'integer':
            metadata.type['convention']['integer'].append(k)
        elif convention['properties'][k]['type'] == 'array':
            metadata.type['convention']['array'].append(k)
        elif convention.get('properties').get(k).get('format') == 'time':
            metadata.type['convention']['number'].append(k)
        if convention.get('properties').get(k).get('ontology'):
            metadata.type['convention']['ontology'].append(k)
    return


def merge_numerics(metadata):
    """Add numeric headers from metadata file to appropriate metadata.fype
    """
    logger.debug('Begin: merge_numerics')
    metadata.type['floats'] = metadata.type['convention']['number'][:]
    for n in metadata.type['numeric_headers']:
        n_not_number = n not in metadata.type['convention']['number']
        n_not_integer = n not in metadata.type['convention']['integer']
        if n_not_number and n_not_integer:
            metadata.type['floats'].append(n)
    return


def list_duplicates(cells):
    """Find duplicates in list for detailed reporting
    """
    # reference https://stackoverflow.com/questions/9835762
    seen = set()
    # save add function to avoid repeated lookups
    seen_add = seen.add
    # list comprehension below chosen for computational efficiency
    # adds all new elements to seen and all other to seen_twice
    # "or" allows seen_add(x) evaluation only if x is not in seen
    seen_twice = set(x for x in cells if x in seen or seen_add(x))
    return list(seen_twice)


def validate_cells_unique(metadata):
    """Check all CellID are unique.

    :return: boolean   True if valid, False otherwise
    """
    valid = False
    if len(metadata.cells) == len(set(metadata.cells)):
        valid = True
    else:
        dups = list_duplicates(metadata.cells)
        metadata.errors['format'].append(
            'Error:  Duplicate CellID(s) in metadata file:'
            f' {", ".join(map(str, dups))}'
        )
    return valid


def collect_ontology_data(data, metadata):
    """Collect unique ontology IDs for ontology validation
    """
    ### function is not yet working, fix extract_convention_types
    logger.debug('Begin: collect_ontology_data')
    local_errors = set()
    for entry in data.keys():
        if entry in metadata.type['convention']['ontology']:
            metadata.ontology[entry]['ontologyID'].add(entry)
            metadata.ontology[entry]['CellID'].add(data['CellID'])
            label_key = entry + '-ontology_label'
        try:
            metadata.ontology[entry]['label'].add(data[label_key])
        except KeyError:
            ### prefer to capture below to errors
            msg = 'Warning: missing key' + label_key
            local_errors.add(msg)
        except UnboundLocalError:
            ### ToDo - figure out best practice for this
            local_errors.add('UnboundLocalError - needs to be addressed')
    return


def process_metadata_row(metadata, convention, line):
    """Read TSV metadata input file row by row

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    :return: row of convention data
    """
    # linter complaining about complexity index, suggestions welcomed
    logger.debug('Begin: process_metadata_input')
    extract_numeric_headers(metadata)
    extract_convention_types(convention, metadata)
    merge_numerics(metadata)
    metadata.headers[0] = 'CellID'
    keys = metadata.headers
    row_info = dict(zip(keys, line))
    for k, v in row_info.items():
        if not v:
            row_info[k] = None
        if k in metadata.type['convention']['integer']:
            try:
                row_info[k] = int(v)
            except ValueError as e:
                print(
                    'Warning: Issue with coercion based on type declaration', e
                )
        elif k in metadata.type['floats']:
            try:
                row_info[k] = float(v)
            except ValueError as e:
                print(
                    'Warning: Issue with coercion based on type declaration', e
                )
        elif k in metadata.type['convention']['array']:
            try:
                row_info[k] = v.split(',')
            except ValueError as e:
                print(
                    'Warning: Issue with coercion based on type declaration', e
                )
    return row_info


def process_metadata_content(metadata, convention):
    """Evaluate TSV metadata input non-ontology errors and ontology info

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    :return: tuple of non-ontology error dict and ontology info dict
            or False if input convention is invalid jsonschema
    """
    logger.debug('Begin: process_metadata_content')
    # this function seems overloaded with its three tasks
    # schema validation, non-ontology errors, ontology info collection
    # the latter two should be done together in the same pass thru the file
    js_errors = defaultdict(list)
    schema = validate_schema(convention)
    if schema:
        line = metadata.extract()
        while line:
            # print('line:', line)
            row = process_metadata_row(metadata, convention, line)
            # print('row:', row)
            metadata.cells.append(row['CellID'])
            collect_ontology_data(row, metadata)
            for error in schema.iter_errors(row):
                js_errors[error.message].append(row['CellID'])
            line = metadata.extract()
        metadata.errors['values'] = js_errors
        validate_cells_unique(metadata)
        return
    else:
        print('Validation failed: Invalid metadata convention')
        return


def print_collected_ontology_data(metadata):
    logger.debug('Begin: print_collected_ontology_data')
    for entry in metadata.ontology.keys():
        print(entry)
        print(
            'Ontology data for entry: (count =',
            len(metadata.ontology[entry]['CellID']), ')'
        )
        print('Uniq ontologyIDs:', metadata.ontology[entry]['ontologyID'])
        print('Uniq labels:', metadata.ontology[entry]['label'])
    return


def report_errors(metadata):
    """Report errors in error object

    :param metadata: cell metadata object
    :return: True if errors are reported, False if no errors to report
    """
    logger.debug('Begin: report_errors')
    errors = False
    for k, v in metadata.errors.items():
        if k == 'format' and v:
            print('Metadata format errors:', v)
            errors = True
        if k == 'values' and v:
            print('Non-ontology metadata errors:')
            for error, cells in v.items():
                print(error, '[ Error count:', len(cells), ']')
                errors = True
    if not errors:
        ### deal with this print statement
        print('No errors detected in input metadata file')
    return errors


# ToDo
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
ensure numeric (TYPE == numeric in TSV file; type number or integer in loom)
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
    with open(args.convention, 'r') as f:
        convention = json.load(f)
    filetsv = args.input_metadata
    metadata = CellMetadata(filetsv)
    print('Validating', filetsv)
    format_valid = metadata.validate_format()
    if format_valid:
        process_metadata_content(metadata, convention)
    report_errors(metadata)
#    print_collected_ontology_data(metadata)
