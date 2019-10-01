"""Validate input metadata TSV file against metadata convention.

DESCRIPTION
This CLI takes a TSV metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using JSON file for Alexandria metadata convention TSV, validate input TSV
$ python3 validate_metadata.py ../../tests/data/AMC_v0.8.json ../../tests/data/metadata_valid.tsv

"""

import argparse
import json
import logging
from collections import defaultdict
import sys
import requests
import urllib.parse as encoder
import re

import colorama
from colorama import Fore
import jsonschema

sys.path.append('..')
try:
    # Used when importing internally and in tests
    from cell_metadata import CellMetadata
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..cell_metadata import CellMetadata

# ToDo set up parameters to adjust log levels
#  logging.basicConfig(level=logging.INFO, filename='app.log', filemode='w',
#   format='%(name)s - %(levelname)s - %(message)s')
logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Ensures normal color for print() output, unless explicitly changed
colorama.init(autoreset=True)


def create_parser():
    """Parse command line values for validate_metadata
    """
    logger.debug('Begin: create_parser')
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # parser.add_argument(
    #     '--output',
    #     '-o',
    #     type=str,
    #     help='Output file name [optional]',
    #     default=None
    # )
    # parser.add_argument(
    #     '--key_id',
    #     '-k'
    #     type=str,
    #     help='Key metadata name for parsing; CellID for metadata, BiosampleID for sample sheets [optional]',
    #     default='CellID'
    # )

    # helper param to create JSON representation of metadata.error
    # as reference output for tests
    parser.add_argument('--issues_json', action='store_true')
    # TODO: make required and modify defaults on the following two parameters after consulting Jon
    parser.add_argument('--file_id', help='MongoDB identifier', default='Mongo_none')
    parser.add_argument(
        '--study_accession', help='SCP study accession', default='SCP_none'
    )
    parser.add_argument('convention', help='Metadata convention JSON file ')
    parser.add_argument('input_metadata', help='Metadata TSV file')
    return parser


def validate_schema(json, metadata):
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
        error_msg = 'Error: Invalid metadata convention'
        metadata.store_validation_issue('error', 'convention', error_msg, repr(e))
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

    CellMetadata.type.convention populated with lists of metadata
    for number, integer, array and ontology
    """
    logger.debug('Begin: extract_convention_types')
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
    # print('convention ontologies', metadata.type['convention']['ontology'])
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
    # 'or' allows seen_add(x) evaluation only if x is not in seen
    seen_twice = set(x for x in cells if x in seen or seen_add(x))
    return list(seen_twice)


def validate_cells_unique(metadata):
    """Check all CellID are unique.

    :return: boolean   True if valid, False otherwise
    """
    valid = False
    # uniq_errors = defaultdict(list)
    if len(metadata.cells) == len(set(metadata.cells)):
        valid = True
    else:
        # TODO value stored incorrectly, will need to fix
        dups = list_duplicates(metadata.cells)
        msg = 'Error:  Duplicate CellID(s) in metadata file'
        # uniq_errors[msg].append(dups)
        # metadata.issues['error']['format'] = uniq_errors
        metadata.issues['error']['format'][msg] = dups
    return valid


def collect_cell_for_ontology(metadatum, data, metadata):
    """Collect ontology info for a single metadatum into CellMetadata.ontology dictionary
    """
    local_errors = set()
    logger.debug('Begin: collect_cell_for_ontology')
    ontology_label = metadatum + '__ontology_label'
    try:
        metadata.ontology[metadatum][(data[metadatum], data[ontology_label])].append(
            data['CellID']
        )
    except KeyError:
        metadata.ontology[metadatum][(data[metadatum])].append(data['CellID'])
    return local_errors


def collect_ontology_data(data, metadata):
    """Collect unique ontology IDs for ontology validation
    """
    logger.debug('Begin: collect_ontology_data')
    for entry in data.keys():
        # print('consider collecting ontology', entry)
        if entry in metadata.type['convention']['ontology']:
            # skip ontologies that are arrays for now
            if entry in metadata.type['convention']['array']:
                pass
                # not collecting metadata for array-based metadata
            else:
                collect_cell_for_ontology(entry, data, metadata)
        else:
            # troubleshooting print statement to visualize skipping metadata that are not of type ontology
            # print(entry, 'not an ontology in', metadata.type['convention']['ontology'])
            pass
    return


def process_metadata_row(metadata, convention, line):
    """Read TSV metadata input file row by row

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    :return: row of convention data
    """
    # linter complaining about complexity index, suggestions welcomed
    logger.debug('Begin: process_metadata_row')
    extract_numeric_headers(metadata)
    extract_convention_types(convention, metadata)
    merge_numerics(metadata)
    metadata.headers[0] = 'CellID'
    keys = metadata.headers
    row_info = dict(zip(keys, line))
    for k, v in row_info.items():
        # explicitly setting empty values to None so missing values for
        # required data are caught at validation
        #        print('processing row metadata:', k)
        if not v:
            row_info[k] = None
        try:
            if k in metadata.type['convention']['integer']:
                row_info[k] = int(v)
            elif k in metadata.type['floats']:
                row_info[k] = float(v)
            elif k in metadata.type['convention']['array']:
                row_info[k] = v.split(',')
        except ValueError:
            error_msg = f'ERROR: value provided does not match convention type declaration for {k}'
            metadata.store_validation_issue(
                'error', 'type', error_msg, [row_info['CellID']]
            )
    return row_info


def collect_jsonschema_errors(metadata, convention):
    """Evaluate TSV metadata input against metadata convention using JSON schema

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    :return: tuple of non-ontology issues dict and ontology info dict
            or False if input convention is invalid JSON schema
    """
    logger.debug('Begin: collect_jsonschema_errors')
    # this function seems overloaded with its three tasks
    # schema validation, non-ontology errors, ontology info collection
    # the latter two should be done together in the same pass thru the file
    js_errors = defaultdict(list)
    schema = validate_schema(convention, metadata)

    if schema:
        metadata.reset_file(2)
        line = metadata.extract()
        row_count = 1
        while line:
            # print('processing row', row_count)
            row_count += 1
            # print('line:', line)
            row = process_metadata_row(metadata, convention, line)
            # print('row:', row)
            metadata.cells.append(row['CellID'])
            collect_ontology_data(row, metadata)
            for error in schema.iter_errors(row):
                js_errors[error.message].append(row['CellID'])
            line = metadata.extract()
        metadata.issues['error']['convention'] = js_errors
        validate_cells_unique(metadata)
        return


def report_issues(metadata):
    """Report issues in CellMetadata.issues dictionary

    :param metadata: cell metadata object
    :return: True if errors are reported, False if no errors to report
    """
    logger.debug('Begin: report_issues')

    has_errors = False
    has_warnings = False
    for issue_type in sorted(metadata.issues.keys()):
        for issue_category, category_dict in metadata.issues[issue_type].items():
            if category_dict:
                print('\n***', issue_category, issue_type, 'list:')
                if issue_type == 'error':
                    color = Fore.RED
                    has_errors = True
                elif issue_type == 'warn':
                    color = Fore.YELLOW
                    has_warnings = True
                for issue_text, cells in category_dict.items():
                    issue_msg = color + issue_text
                    if cells:
                        print(f'{issue_msg} [ Error count: {len(cells)} ]')
                    else:
                        print(issue_msg)
    if not has_errors and not has_warnings:
        # deal with this print statement
        print('No errors or warnings detected for input metadata file')
    return has_errors


def exit_if_errors(metadata):
    """Determine if CellMetadata.issues has errors

    :param metadata: cell metadata object
    :return: Exit with error code 1 if errors are reported, False if no errors
    """
    logger.debug('Begin: exit_if_errors')

    errors = False
    for error_type in metadata.issues.keys():
        for error_category, category_dict in metadata.issues[error_type].items():
            if category_dict:
                if error_type == 'error':
                    errors = True
    if errors:
        exit(1)
    return errors


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


def retrieve_ontology_term(convention_url, ontology_id):
    """Retrieve an individual term from an ontology
    :param ontology_term: term to query for in matching ontology
    :return: JSON payload of ontology of ontology term, or None
    """
    OLS_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/"
    convention_ontology = retrieve_ontology(convention_url)
    # separate ontology shortname from term ID number
    # valid separators are underscore and colon (used by HCA)
    try:
        ontology_shortname, term_id = re.split('[_:]', ontology_id)
    # when ontolgyID is malformed and has no separator -> ValueError
    # when ontologyID value is empty string -> TypeError
    except (ValueError, TypeError) as error:
        print("Exception:", error)
        print('Problem with provided ontology ID', ontology_id)
        return None
    metadata_url = OLS_BASE_URL + ontology_shortname
    metadata_ontology = retrieve_ontology(metadata_url)
    if convention_ontology and metadata_ontology:
        base_term_uri = metadata_ontology['config']['baseUris'][0]
        query_iri = encode_term_iri(term_id, base_term_uri)
        term_url = convention_ontology['_links']['terms']['href'] + '/' + query_iri
        # print(term_url)
        response = requests.get(term_url)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    elif not metadata_ontology:
        print(
            f'No result from EBI OLS for provided ontology shortname \"{ontology_shortname}\"'
        )
    else:
        print(f'encountered issue retrieving {convention_url} or {ontology_shortname}')
        return None


def encode_term_iri(term_id, base_uri):
    """Double url-encode a term Internationalized Resource Identifier (IRI) for querying OLS ontologies

    :param term: ontology term
    :param base_uri: base term URI for corresponding ontology
    :return: double url-encoded ontology term IRI
    """
    query_uri = base_uri + term_id
    encoded_iri = encoder.quote_plus(encoder.quote_plus(query_uri))
    return encoded_iri


def validate_collected_ontology_data(metadata, convention):
    logger.debug('Begin: validate_collected_ontology_data')
    for entry in metadata.ontology.keys():
        # print(entry)
        ontology_url = convention['properties'][entry]['ontology']
        try:
            for ontology_id, ontology_label in metadata.ontology[entry].keys():
                # print(f'checking {ontology_id} and {ontology_label}')
                matching_term = retrieve_ontology_term(ontology_url, ontology_id)
                if matching_term:
                    # print("Found matching for " + ontology_id)
                    if matching_term['label'] != ontology_label:
                        # print(
                        #     f'checking if {ontology_label} matches: {matching_term["label"]}'
                        # )
                        try:
                            error_msg = (
                                f'ontology_label \"{ontology_label}\" '
                                f'does not match \"{matching_term["label"]}\",'
                            )
                            metadata.store_validation_issue(
                                'error',
                                'ontology',
                                error_msg,
                                metadata.ontology[entry][(ontology_id, ontology_label)],
                            )
                        except TypeError:
                            error_msg = f'No description found for {ontology_id} in ontology lookup'
                            metadata.store_validation_issue(
                                'error',
                                'ontology',
                                error_msg,
                                metadata.ontology[entry][(ontology_id, ontology_label)],
                            )
                    # else:
                    #     print("Ontology label matches: " + matching_term['label'])
                else:
                    if ontology_id:
                        error_msg = f'No match found for {ontology_id}'
                        metadata.store_validation_issue(
                            'error',
                            'ontology',
                            error_msg,
                            metadata.ontology[entry][(ontology_id, ontology_label)],
                        )
                    else:
                        error_msg = f'No ontology_id provided for {ontology_label}'
                        print(error_msg)
                        # FIX issue with adding this error to issues, causes
                        # "RuntimeError: dictionary changed size during iteration"
                        # referencing "for ontology_id, ontology_label in metadata.ontology[entry].keys():"
                        # metadata.store_validation_issue(
                        #     'error',
                        #     'ontology',
                        #     error_msg,
                        #     metadata.ontology[entry][(ontology_label)],
                        # )
        except ValueError:
            for ontology_id in metadata.ontology[entry].keys():
                matching_term = retrieve_ontology_term(ontology_url, ontology_id)
                if not matching_term:
                    error_msg = f'No match found for {ontology_id}'
                    metadata.store_validation_issue(
                        'error',
                        'ontology',
                        error_msg,
                        metadata.ontology[entry][(ontology_id)],
                    )
                # else:
                #     print("Found matching for " + ontology_id)
                error_msg = (
                    f'Warning: no ontology label supplied in metadata file for '
                    f'\"{ontology_id}\" - no cross-check for data entry error possible'
                )
                metadata.store_validation_issue('warn', 'ontology', error_msg)
    return


def serialize_issues(metadata):
    with open('issues.json', 'w') as jsonfile:
        json.dump(metadata.issues, jsonfile)


# ToDo
"""
WAIT: handle array data types

"""
"""
DEFER: Things to check before pass intended data types to Firestore
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
    metadata = CellMetadata(filetsv, args.file_id, args.study_accession)
    print('Validating', filetsv)
    format_valid = metadata.validate_format()
    collect_jsonschema_errors(metadata, convention)
    validate_collected_ontology_data(metadata, convention)
    if args.issues_json:
        serialize_issues(metadata)
    report_issues(metadata)
    exit_if_errors(metadata)
