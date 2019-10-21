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

    # helper param to create JSON representation of metadata.issues
    # to generate reference output for tests
    parser.add_argument('--issues_json', action='store_true')
    # validate_metadata.py CLI only for dev, bogus defaults below shouldn't propagate
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
    except jsonschema.SchemaError:
        error_msg = 'Invalid metadata convention file, cannot validate metadata.'
        metadata.store_validation_issue('error', 'convention', error_msg)
        return None


def is_array_metadata(convention, metadatum):
    """Check if metadata is array type from metadata convention

    :param convention: dict representation of metadata convention
    :param metadatum: name of metadatum

    :return:
    """
    logger.debug('Begin: is_array_metadata')
    try:
        type_lookup = convention['properties'][metadatum]['type']
        if type_lookup == 'array':
            return True
        else:
            return False
    except KeyError:
        return False


def is_ontology_metadata(convention, metadatum):
    """Check if metadata is ontology from metadata convention

    :param convention: dict representation of metadata convention
    :param metadatum: name of metadatum

    :return:
    """
    logger.debug('Begin: is_ontology_metadata')
    try:
        if convention['properties'][metadatum]['ontology']:
            return True
        else:
            return False
    except KeyError:
        return False


def lookup_metadata_type(convention, metadatum):
    """Look up metadata type from metadata convention

    :param convention: dict representation of metadata convention
    :param metadatum: name of metadatum

    :return:
    """
    logger.debug('Begin: lookup_metadata_type')
    try:
        if is_array_metadata(convention, metadatum):
            type_lookup = convention['properties'][metadatum]['items']['type']
        else:
            type_lookup = convention['properties'][metadatum]['type']
        return type_lookup
    except KeyError:
        return None


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
    if len(metadata.cells) == len(set(metadata.cells)):
        valid = True
    else:
        dups = list_duplicates(metadata.cells)
        error_msg = 'Duplicate CellID(s) in metadata file'
        metadata.store_validation_issue('error', 'format', error_msg, dups)
    return valid


def collect_cell_for_ontology(metadatum, row_data, metadata, array=False):
    """Collect ontology info for a single metadatum into CellMetadata.ontology dictionary
    """
    logger.debug('Begin: collect_cell_for_ontology')
    if metadatum.endswith('__unit'):
        ontology_label = metadatum + '_label'
    else:
        ontology_label = metadatum + '__ontology_label'
    if array:
        try:
            ontology_dict = dict(zip(row_data[metadatum], row_data[ontology_label]))
            for id, label in ontology_dict.items():
                metadata.ontology[metadatum][(id, label)].append(row_data['CellID'])
        except TypeError:
            for id in row_data[metadatum]:
                metadata.ontology[metadatum][(id)].append(row_data['CellID'])
    else:
        try:
            metadata.ontology[metadatum][
                (row_data[metadatum], row_data[ontology_label])
            ].append(row_data['CellID'])
        except KeyError:
            metadata.ontology[metadatum][(row_data[metadatum])].append(
                row_data['CellID']
            )
    return


def collect_ontology_data(row_data, metadata, convention):
    """Collect unique ontology IDs for ontology validation
    """
    logger.debug('Begin: collect_ontology_data')
    for entry in row_data.keys():
        if is_ontology_metadata(convention, entry):
            if is_array_metadata(convention, entry):
                collect_cell_for_ontology(entry, row_data, metadata, array=True)
            else:
                collect_cell_for_ontology(entry, row_data, metadata)
    return


def cast_boolean_type(value):
    if value.lower() == 'true':
        return True
    elif value.lower() == 'false':
        return False
    else:
        return value


def cast_integer_type(value):
    return int(value)


def cast_float_type(value):
    return float(value)


def return_without_cast(value):
    return value


def cast_metadata_type(metadatum, value, row_info, convention, metadata):
    metadata_types = {
        'number': cast_float_type,
        'boolean': cast_boolean_type,
        'integer': cast_integer_type,
        'string': return_without_cast,
    }
    if is_array_metadata(convention, metadatum):
        cast_values = []
        try:
            # splitting on pipe character for array data, valid for Sarah's
            # programmatically generated SCP TSV metadata files. When ingesting
            # files that support array-based metadata navtively (eg. loom,
            # anndata etc) splitting on pipe may become problematic
            for element in value.split('|'):
                cast_element = metadata_types.get(
                    lookup_metadata_type(convention, metadatum)
                )(element)
                cast_values.append(cast_element)
            row_info[metadatum] = cast_values
        except ValueError:
            error_msg = (
                f'{metadatum}: "{element}" in "{value}" does not match expected type'
            )
            metadata.store_validation_issue(
                'error', 'type', error_msg, [row_info['CellID']]
            )
        # This exception should only trigger if a single-value boolean array
        # metadata is being cast - the value needs to be passed as an array,
        # it is already boolean via Pandas' inference processes
        except AttributeError:
            row_info[metadatum] = [value]
    else:
        try:
            cast_value = metadata_types.get(
                lookup_metadata_type(convention, metadatum)
            )(value)
            row_info[metadatum] = cast_value
        except ValueError:
            error_msg = f'{metadatum}: "{value}" does not match expected type'
            metadata.store_validation_issue(
                'error', 'type', error_msg, [row_info['CellID']]
            )
    return row_info


def compare_type_annots_to_convention(metadata, convention):
    """Check metadata type annotation consistent with metadata convention type

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    """
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    type_annots = metadata.file.columns.get_level_values(1).tolist()
    metadata_names[0] = 'CellID'
    type_annots[0] = 'group'
    metadata_annots = dict(zip(metadata_names, type_annots))
    annot_equivalents = {
        'numeric': ['number', 'integer'],
        'group': ['boolean', 'string'],
    }
    for metadatum, annot in metadata_annots.items():
        convention_type = lookup_metadata_type(convention, metadatum)
        try:
            if convention_type and convention_type not in annot_equivalents.get(annot):
                for k, v in annot_equivalents.items():
                    if convention_type in v:
                        expected = k
                error_msg = (
                    f'{metadatum}: "{annot}" annotation in metadata file conflicts with metadata convention. '
                    f'Convention expects "{expected}" values.'
                )
                metadata.store_validation_issue('error', 'type', error_msg)
        except TypeError:
            for k, v in annot_equivalents.items():
                if convention_type in v:
                    expected = k
            error_msg = (
                f'{metadatum}: "{annot}" annotation in metadata file disagrees with metadata convention. '
                f'Convention expects "{expected}" annotation.'
            )
            metadata.store_validation_issue('error', 'type', error_msg)


def process_metadata_row(metadata, convention, line):
    """Read TSV metadata input file row by row

    :param metadata: cell metadata object
    :param convention: dict representation of metadata convention
    :return: row of convention data
    """
    logger.debug('Begin: process_metadata_row')
    # extract first row of metadata file from pandas array as python list
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    metadata_names[0] = 'CellID'
    row_info = dict(zip(metadata_names, line))
    for k, v in row_info.items():
        row_info = cast_metadata_type(k, v, row_info, convention, metadata)
    return row_info


def collect_jsonschema_errors(metadata, convention):
    """Evaluate metadata input against metadata convention using JSON schema

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
        compare_type_annots_to_convention(metadata, convention)
        rows = metadata.yield_by_row()
        line = next(rows)
        while line:
            row = process_metadata_row(metadata, convention, line)
            metadata.cells.append(row['CellID'])
            collect_ontology_data(row, metadata, convention)
            for error in schema.iter_errors(row):
                try:
                    error.message = error.path[0] + ': ' + error.message
                except IndexError:
                    pass
                js_errors[error.message].append(row['CellID'])
            try:
                line = next(rows)
            except StopIteration:
                break
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
    OLS_BASE_URL = 'https://www.ebi.ac.uk/ols/api/ontologies/'
    convention_ontology = retrieve_ontology(convention_url)
    # separate ontology shortname from term ID number
    # valid separators are underscore and colon (used by HCA)
    try:
        ontology_shortname, term_id = re.split('[_:]', ontology_id)
    # when ontolgyID is malformed and has no separator -> ValueError
    # when ontologyID value is empty string -> TypeError
    except (ValueError, TypeError):
        return None
    metadata_url = OLS_BASE_URL + ontology_shortname
    metadata_ontology = retrieve_ontology(metadata_url)
    if convention_ontology and metadata_ontology:
        base_term_uri = metadata_ontology['config']['baseUris'][0]
        query_iri = encode_term_iri(term_id, base_term_uri)
        term_url = convention_ontology['_links']['terms']['href'] + '/' + query_iri
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
    """Evaluate collected ontology_id, ontology_label info in
    CellMetadata.ontology dictionary by querying EBI OLS for
    validity of ontology_id and cross-check that input ontology_label
    matches ontology_label from EBI OLS lookup
    """
    logger.debug('Begin: validate_collected_ontology_data')
    for entry in metadata.ontology.keys():
        ontology_url = convention['properties'][entry]['ontology']
        try:
            for ontology_id, ontology_label in metadata.ontology[entry].keys():
                matching_term = retrieve_ontology_term(ontology_url, ontology_id)
                if matching_term:
                    if matching_term['label'] != ontology_label:
                        try:
                            error_msg = (
                                f'{entry}: input ontology_label \"{ontology_label}\" '
                                f'does not match EBI OLS lookup \"{matching_term["label"]}\",'
                            )
                            metadata.store_validation_issue(
                                'error',
                                'ontology',
                                error_msg,
                                metadata.ontology[entry][(ontology_id, ontology_label)],
                            )
                        except TypeError:
                            error_msg = f'{entry}: No description found for {ontology_id} in ontology lookup'
                            metadata.store_validation_issue(
                                'error',
                                'ontology',
                                error_msg,
                                metadata.ontology[entry][(ontology_id, ontology_label)],
                            )
                # handle case where EBI OLS has no match result
                else:
                    # empty cells for ontology_id and ontology_label now nan when using pandas for ingest
                    if ontology_id:
                        error_msg = f'{entry}: No match found for {ontology_id}'
                        metadata.store_validation_issue(
                            'error',
                            'ontology',
                            error_msg,
                            metadata.ontology[entry][(ontology_id, ontology_label)],
                        )
                    # conditions for this else clause should no longer occur (with pandas)
                    # if this clause is triggered, a RuntimeError should result
                    else:
                        error_msg = (
                            f'{entry}: No ontology_id provided for {ontology_label}'
                        )
                        metadata.store_validation_issue(
                            'error',
                            'ontology',
                            error_msg,
                            metadata.ontology[entry][(ontology_label)],
                        )
        # handle case where no ontology_label provided
        except (TypeError, ValueError):
            for ontology_id in metadata.ontology[entry].keys():
                matching_term = retrieve_ontology_term(ontology_url, ontology_id)
                if not matching_term:
                    error_msg = f'{entry}: No match found for {ontology_id}'
                    metadata.store_validation_issue(
                        'error',
                        'ontology',
                        error_msg,
                        metadata.ontology[entry][(ontology_id)],
                    )
                error_msg = (
                    f'{entry}: no ontology label supplied in metadata file for '
                    f'\"{ontology_id}\" - no cross-check for data entry error possible'
                )
                metadata.store_validation_issue('warn', 'ontology', error_msg)
    return


def confirm_uniform_units(metadata, convention):
    """Check that any unit metadata are uniform within study
    Note: refactoring may be needed if metadata files are chunked
    """
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    for name in metadata_names:
        if name.endswith('__unit'):
            if metadata.file[name].nunique(dropna=False).values[0] != 1:
                error_msg = (
                    f'{name}: values for each unit metadata required to be uniform'
                )
                metadata.store_validation_issue('error', 'convention', error_msg)


def serialize_issues(metadata):
    """Write collected issues to json file
    """
    with open('issues.json', 'w') as jsonfile:
        json.dump(metadata.issues, jsonfile, indent=2)


def validate_input_metadata(metadata, convention):
    """Wrapper function to run validation functions
    """
    collect_jsonschema_errors(metadata, convention)
    validate_collected_ontology_data(metadata, convention)
    confirm_uniform_units(metadata, convention)


if __name__ == '__main__':
    args = create_parser().parse_args()
    arguments = vars(args)
    with open(args.convention, 'r') as f:
        convention = json.load(f)
    filetsv = args.input_metadata
    metadata = CellMetadata(
        filetsv, args.file_id, args.study_accession, open_as='dataframe'
    )
    print('Validating', filetsv)
    validate_input_metadata(metadata, convention)
    if args.issues_json:
        serialize_issues(metadata)
    report_issues(metadata)
    exit_if_errors(metadata)
