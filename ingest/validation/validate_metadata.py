"""Validate input metadata TSV file against metadata convention.

DESCRIPTION
This CLI takes a TSV metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using JSON file for Alexandria metadata convention TSV, validate input TSV
$ python3 validate_metadata.py ../../tests/data/AMC_v1.1.3.json ../../tests/data/valid_no_array_v1.1.3.tsv \
          --study-accession SCP123 --study-id 5dfa6718421aa90fea085476 --study-file-id 5e27451e2209c211b1e7c9cc

"""

import argparse
import json
import logging
from collections import defaultdict
import sys
import requests
import urllib.parse as encoder
import re
import os
import numbers
import time
import backoff

import colorama
from colorama import Fore
import jsonschema
from google.cloud import bigquery

sys.path.append('..')
try:
    # Used when importing internally and in tests
    from cell_metadata import CellMetadata
    from monitor import setup_logger
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..cell_metadata import CellMetadata
    from .monitor import setup_logger


info_logger = setup_logger(__name__, "info.txt")

error_logger = setup_logger(
    __name__ + "_errors", "errors.txt", level=logging.ERROR, format=''
)

# Ensures normal color for print() output, unless explicitly changed
colorama.init(autoreset=True)

# Configure maximum number of seconds to spend & total attempts at external HTTP requests to services, e.g. OLS
MAX_HTTP_REQUEST_TIME = 120
MAX_HTTP_ATTEMPTS = 8


def create_parser():
    """Parse command line values for validate_metadata
    """
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
    parser.add_argument('--issues-json', action='store_true')
    # helper param to create JSON representation of convention metadata
    # to generate json for bigquery testing
    parser.add_argument('--bq-json', action='store_true')
    # overwrite existing output
    parser.add_argument('--force', action='store_true')
    # test BigQuery upload functions
    parser.add_argument('--upload', action='store_true')
    # validate_metadata.py CLI only for dev, bogus defaults below shouldn't propagate
    parser.add_argument(
        '--study-id',
        help='MongoDB study identifier',
        default='6d276a50421aa9117c982846',
    )
    parser.add_argument(
        '--study-file-id', help='MongoDB file identifier', default='abc123'
    )
    parser.add_argument(
        '--study-accession', help='SCP study accession', default='SCP888'
    )
    parser.add_argument(
        '--bq-dataset', help='BigQuery dataset identifier', default='cell_metadata'
    )
    parser.add_argument(
        '--bq-table', help='BigQuery table identifier', default='alexandria_convention'
    )
    parser.add_argument('convention', help='Metadata convention JSON file ')
    parser.add_argument('input_metadata', help='Metadata TSV file')
    return parser


def validate_schema(json, metadata):
    """Check validity of metadata convention as JSON schema.
    if valid, return jsonschema validator object else return None
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
    """
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
    """
    try:
        return bool(convention['properties'][metadatum]['ontology'])
    except KeyError:
        return False


def lookup_metadata_type(convention, metadatum):
    """Look up metadata type from metadata convention
    """
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
    for entry in row_data.keys():
        if is_ontology_metadata(convention, entry):
            if is_array_metadata(convention, entry):
                collect_cell_for_ontology(entry, row_data, metadata, array=True)
            else:
                collect_cell_for_ontology(entry, row_data, metadata)
    return


def compare_type_annots_to_convention(metadata, convention):
    """Check if metadata type annotation is consistent with metadata convention type
    """
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    type_annots = metadata.file.columns.get_level_values(1).tolist()
    # if input was TSV metadata file, SCP format requires 'NAME' for the first
    # column which is expected to be CellID, primarily based on loom convention
    # would do conditional (eg. if metadata_names[0].upper() == 'NAME':)
    # but subsequent code expects 'CellID' even if header check fails and first
    # word in header is not some form of 'NAME' and script breaks
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
            if '.' in annot:
                # duplicated metadata header name identified in validate_unique_header
                # detection of invalid type annotation is side effect of Pandas
                pass
            elif 'Unnamed' in annot:
                # missing type annotation also detected in validate_against_header_count
                # invalid type annotation is side effect of Pandas
                error_msg = (
                    f'{metadatum}: missing TYPE annotation in metadata file. '
                    f'Convention expects "{expected}" annotation.'
                )
                metadata.store_validation_issue('error', 'format', error_msg)
            else:
                error_msg = (
                    f'{metadatum}: invalid "{annot}" annotation in metadata file. '
                    f'Convention expects "{expected}" annotation.'
                )
                metadata.store_validation_issue('error', 'type', error_msg)


def cast_boolean_type(value):
    """Cast metadata value as boolean, if castable
    """
    if isinstance(value, bool):
        return value
    elif str(value).lower() == 'true':
        return True
    elif str(value).lower() == 'false':
        return False
    else:
        raise ValueError(f'cannot cast {value} as boolean')


def cast_integer_type(value):
    """Cast metadata value as integer
    """
    return int(value)


def cast_float_type(value):
    """Cast metadata value as float
    """
    return float(value)


def cast_string_type(value):
    """Cast string type per convention where Pandas autodetected a number
    """
    if isinstance(value, numbers.Number):
        return str(value)
    else:
        return value

def regularize_ontologyID(value):
    """Regularize ontologyIDs for storage with underscore format
    """
    try:
        ontology_shortname, term_id = re.split('[_:]', value)
        value = ontology_shortname + '_' + term_id
        return value
    # when ontolgyID is malformed and has no separator -> ValueError
    # when ontologyID value is empty string -> TypeError
    except (ValueError, TypeError):
        return value

def cast_metadata_type(metadatum, value, id_for_error_detail, convention, metadata):
    """for metadatum, lookup expected type by metadata convention
        and cast value as appropriate type for validation
    """
    cast_metadata = {}
    metadata_types = {
        'number': cast_float_type,
        'boolean': cast_boolean_type,
        'integer': cast_integer_type,
        'string': cast_string_type,
    }
    if is_array_metadata(convention, metadatum):
        cast_values = []
        try:
            # splitting on pipe character for array data, valid for Sarah's
            # programmatically generated SCP TSV metadata files. When ingesting
            # files that support array-based metadata navtively (eg. loom,
            # anndata etc) splitting on pipe may become problematic
            for element in value.split('|'):
                if 'ontology' in convention['properties'][metadatum]:
                    element = regularize_ontologyID(element)
                cast_element = metadata_types.get(
                    lookup_metadata_type(convention, metadatum)
                )(element)
                cast_values.append(cast_element)
            cast_metadata[metadatum] = cast_values
        except ValueError:
            error_msg = (
                f"{metadatum}: '{element}' in '{value}' does not match "
                f"expected '{lookup_metadata_type(convention, metadatum)}' type"
            )
            metadata.store_validation_issue(
                'error', 'type', error_msg, [id_for_error_detail]
            )
        # This exception should only trigger if a single-value boolean array
        # metadata is being cast - the value needs to be passed as an array,
        # it is already boolean via Pandas' inference processes
        except AttributeError:
            if 'ontology' in convention['properties'][metadatum]:
                value = regularize_ontologyID(value)
            cast_metadata[metadatum] = [metadata_types.get(
                lookup_metadata_type(convention, metadatum)
            )(value)]
    else:
        if 'ontology' in convention['properties'][metadatum]:
            value = regularize_ontologyID(value)
        try:
            cast_value = metadata_types.get(
                lookup_metadata_type(convention, metadatum)
            )(value)
            cast_metadata[metadatum] = cast_value
        except ValueError:
            error_msg = f'{metadatum}: "{value}" does not match expected type'
            metadata.store_validation_issue(
                'error', 'type', error_msg, [id_for_error_detail]
            )
        # particular metadatum is not in convention, metadata does not need
        # to be added to new_row for validation, return empty dictionary
        except TypeError:
            return {}
    return cast_metadata


def process_metadata_row(metadata, convention, line):
    """Process metadata row by row
    returns processed row of convention data as dict
    """
    # extract first row of metadata from pandas array as python list
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    # if input was TSV metadata file, SCP format requires 'NAME' for the first
    # column which is expected to be CellID, primarily based on loom convention
    # would do conditional (eg. if metadata_names[0].upper() == 'NAME':)
    # but subsequent code expects 'CellID' even if header check fails and first
    # word in header is not some form of 'NAME' and script breaks
    metadata_names[0] = 'CellID'
    row_info = dict(zip(metadata_names, line))
    processed_row = {}
    for k, v in row_info.items():
        processed_row.update(
            cast_metadata_type(k, v, row_info['CellID'], convention, metadata)
        )
    return processed_row


def collect_jsonschema_errors(metadata, convention, bq_json=None):
    """Evaluate metadata input against metadata convention using JSON schema
    returns False if input convention is invalid JSON schema
    """
    # this function seems overloaded with its three tasks
    # schema validation, non-ontology errors, ontology info collection
    # the latter two should be done together in the same pass thru the file
    js_errors = defaultdict(list)
    schema = validate_schema(convention, metadata)

    if bq_json:
        bq_filename = str(metadata.study_file_id) + '.json'
        # truncate JSON file so data from serialize_bq starts with an empty file
        fh = open(bq_filename, 'w')
        fh.close()
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
            if bq_json:
                # add non-convention, SCP-required, metadata for BigQuery
                row['study_accession'] = metadata.study_accession
                row['file_id'] = str(metadata.study_file_id)
                serialize_bq(row, bq_filename)
            try:
                line = next(rows)
            except StopIteration:
                break
        metadata.issues['error']['convention'] = js_errors
        validate_cells_unique(metadata)
    else:
        return False


def record_issue(errfile, warnfile, issue_type, msg):
    """print issue to console with coloring and
    writes issues to appropriate issue file
    """

    if issue_type == 'error':
        errfile.write(msg + '\n')
        error_logger.error(msg)
        color = Fore.RED
    elif issue_type == 'warn':
        warnfile.write(msg + '\n')
        color = Fore.YELLOW
    else:
        color = ''
    console_msg = color + msg
    print(console_msg)


def report_issues(metadata):
    """Report issues in CellMetadata.issues dictionary
    returns True if errors are reported, False if no errors to report
    """

    info_logger.info(
        f'Checking for validation issues',
        extra={'study_id': metadata.study_id, 'duration': None},
    )

    error_file = open('scp_validation_errors.txt', 'w')
    warn_file = open('scp_validation_warnings.txt', 'w')
    has_errors = False
    has_warnings = False
    for issue_type in sorted(metadata.issues.keys()):
        for issue_category, category_dict in metadata.issues[issue_type].items():
            if category_dict:
                category_header = f'\n*** {issue_category} {issue_type} list:'
                record_issue(error_file, warn_file, issue_type, category_header)
                if issue_type == 'error':
                    has_errors = True
                elif issue_type == 'warn':
                    has_warnings = True
                for issue_text, cells in category_dict.items():
                    if cells:
                        issue_msg = f'{issue_text} [ Error count: {len(cells)} ]'
                        record_issue(error_file, warn_file, issue_type, issue_msg)
                    else:
                        record_issue(error_file, warn_file, issue_type, issue_text)
    if not has_errors and not has_warnings:
        no_issues = 'No errors or warnings detected for input metadata file'
        record_issue(error_file, warn_file, None, no_issues)
    error_file.close()
    warn_file.close()
    return has_errors


def exit_if_errors(metadata):
    """Determine if CellMetadata.issues has errors
    Exit with error code 1 if errors are reported, return False if no errors
    """
    errors = False
    for error_type in metadata.issues.keys():
        for error_category, category_dict in metadata.issues[error_type].items():
            if category_dict:
                if error_type == 'error':
                    errors = True
    if errors:
        exit(1)
    return errors


def backoff_handler(details):
    """Handler function to log backoff attempts when querying OLS"""
    info_logger.debug(
        "Backing off {wait:0.1f} seconds after {tries} tries "
        "calling function {target} with args {args} and kwargs "
        "{kwargs}".format(**details),
        extra={'study_id': None, 'duration': None},
    )


# Attach exponential backoff to external HTTP requests
@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=MAX_HTTP_REQUEST_TIME,
    max_tries=MAX_HTTP_ATTEMPTS,
    on_backoff=backoff_handler,
    logger="info_logger",
)
def retrieve_ontology(ontology_url):
    """Retrieve an ontology listing from EBI OLS
    returns JSON payload of ontology, or None if unsuccessful
    """
    # add timeout to prevent request from hanging indefinitely
    response = requests.get(ontology_url, timeout=60)
    # inserting sleep to minimize 'Connection timed out' error with too many concurrent requests
    time.sleep(0.25)
    if response.status_code == 200:
        return response.json()
    else:
        return None


# Attach exponential backoff to external HTTP requests
@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=MAX_HTTP_REQUEST_TIME,
    max_tries=MAX_HTTP_ATTEMPTS,
    on_backoff=backoff_handler,
    logger="info_logger",
)
def retrieve_ontology_term(convention_url, ontology_id, ontologies):
    """Retrieve an individual term from an ontology
    returns JSON payload of ontology, or None if unsuccessful
    Will store any retrieved ontologies for faster validation of downstream terms
    """
    OLS_BASE_URL = 'https://www.ebi.ac.uk/ols/api/ontologies/'
    # separate ontology shortname from term ID number
    # valid separators are underscore and colon (used by HCA)
    try:
        ontology_shortname, term_id = re.split('[_:]', ontology_id)
    # when ontolgyID is malformed and has no separator -> ValueError
    # when ontologyID value is empty string -> TypeError
    except (ValueError, TypeError):
        return None
    # check if we have already retrieved this ontology reference
    if ontology_shortname not in ontologies:
        metadata_url = OLS_BASE_URL + ontology_shortname
        ontologies[ontology_shortname] = retrieve_ontology(metadata_url)
    metadata_ontology = ontologies[ontology_shortname]
    # check if the ontology parsed from the term is the same ontology defined in the convention
    # if so, skip the extra call to OLS; otherwise, retrieve the convention-defined ontology for term lookups
    if metadata_ontology is not None:
        reference_url = metadata_ontology['_links']['self']['href']
        if reference_url.lower() == convention_url.lower():
            convention_ontology = metadata_ontology.copy()
        else:
            convention_shortname = extract_terminal_pathname(convention_url)
            convention_ontology = retrieve_ontology(convention_url)
            ontologies[convention_shortname] = convention_ontology
    else:
        convention_ontology = (
            None
        )  # we did not get a metadata_ontology, so abort the check
    if convention_ontology and metadata_ontology:
        base_term_uri = metadata_ontology['config']['baseUris'][0]
        query_iri = encode_term_iri(term_id, base_term_uri)
        term_url = convention_ontology['_links']['terms']['href'] + '/' + query_iri
        # add timeout to prevent request from hanging indefinitely
        response = requests.get(term_url, timeout=60)
        # inserting sleep to minimize 'Connection timed out' error with too many concurrent requests
        time.sleep(0.25)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    elif not metadata_ontology:
        error_msg = f'No result from EBI OLS for provided ontology shortname \"{ontology_shortname}\"'
        print(error_msg)
        info_logger.info(error_msg, extra={'study_id': None, 'duration': None})
    else:
        error_msg = (
            f'encountered issue retrieving {convention_url} or {ontology_shortname}'
        )
        print(error_msg)
        info_logger.info(error_msg, extra={'study_id': None, 'duration': None})
        return None


def encode_term_iri(term_id, base_uri):
    """Double url-encode a term Internationalized Resource Identifier (IRI) for querying OLS ontologies
    returns double url-encoded ontology term IRI
    """
    query_uri = base_uri + term_id
    encoded_iri = encoder.quote_plus(encoder.quote_plus(query_uri))
    return encoded_iri


def extract_terminal_pathname(url):
    """Extract the last path segment from a URL
    """
    return list(filter(None, url.split("/"))).pop()


def validate_collected_ontology_data(metadata, convention):
    """Evaluate collected ontology_id, ontology_label info in
    CellMetadata.ontology dictionary by querying EBI OLS for
    validity of ontology_id and cross-check that input ontology_label
    matches ontology_label from EBI OLS lookup
    """
    # container to store references to retrieved ontologies for faster validation
    stored_ontologies = {}
    for entry in metadata.ontology.keys():
        ontology_url = convention['properties'][entry]['ontology']
        try:
            for ontology_id, ontology_label in metadata.ontology[entry].keys():
                matching_term = retrieve_ontology_term(
                    ontology_url, ontology_id, stored_ontologies
                )
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
                        error_msg = f'{entry}: No match found in EBI OLS for provided ontology ID: {ontology_id}'
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
                matching_term = retrieve_ontology_term(
                    ontology_url, ontology_id, stored_ontologies
                )
                if not matching_term:
                    error_msg = f'{entry}: No match found in EBI OLS for provided ontology ID: {ontology_id}'
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
        except requests.exceptions.RequestException as err:
            error_msg = f'External service outage connecting to {ontology_url} when querying {ontology_id}:{ontology_label}: {err}'
            error_logger.error(error_msg)
            metadata.store_validation_issue('error', 'ontology', error_msg)
            # immediately return as validation cannot continue
            return
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


def serialize_bq(bq_dict, filename='bq.json'):
    """Write metadata collected for validation to json file
    BigQuery requires newline delimited json objects
    """
    data = json.dumps(bq_dict)
    with open(filename, 'a') as jsonfile:
        jsonfile.write(data + '\n')


def serialize_issues(metadata):
    """Write collected issues to json file
    """
    with open('issues.json', 'w') as jsonfile:
        json.dump(metadata.issues, jsonfile, indent=2)


def review_metadata_names(metadata):
    """Check metadata names for disallowed characters
    """
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    for name in metadata_names:
        allowed_char = re.compile('^[A-Za-z0-9_]+$')
        if not allowed_char.match(name):
            error_msg = (
                f'{name}: only alphanumeric characters and underscore '
                f'allowed in metadata name'
            )
            metadata.store_validation_issue('error', 'metadata_name', error_msg)


def validate_input_metadata(metadata, convention, bq_json=None):
    """Wrapper function to run validation functions
    """
    collect_jsonschema_errors(metadata, convention, bq_json)
    review_metadata_names(metadata)
    validate_collected_ontology_data(metadata, convention)
    confirm_uniform_units(metadata, convention)


def push_metadata_to_bq(metadata, ndjson, dataset, table):
    """upload local NDJSON to BigQuery
    """
    client = bigquery.Client()
    dataset_ref = client.dataset(dataset)
    table_ref = dataset_ref.table(table)
    job_config = bigquery.LoadJobConfig()
    job_config.write_disposition = bigquery.WriteDisposition.WRITE_APPEND
    job_config.source_format = bigquery.SourceFormat.NEWLINE_DELIMITED_JSON
    try:
        with open(ndjson, 'rb') as source_file:
            job = client.load_table_from_file(
                source_file, table_ref, job_config=job_config
            )
        job.result()  # Waits for table load to complete.
        info_msg = f'Metadata uploaded to BigQuery. ({job.output_rows} rows)'
        print(info_msg)
        info_logger.info(
            info_msg, extra={'study_id': metadata.study_id, 'duration': None}
        )
    # Unable to intentionally trigger a failed BigQuery upload
    # please add print statement below to error logging so
    # error handling can be updated when better understood
    except Exception as e:
        print(e)
        info_logger.info(e, extra={'study_id': metadata.study_id, 'duration': None})
        return 1
    if job.output_rows != len(metadata.cells):
        error_msg = f'BigQuery upload error: upload ({job.output_rows} rows) does not match number of cells in file, {len(metadata.cells)} cells'
        print(error_msg)
        error_logger.error(error_msg)
        return 1
    os.remove(ndjson)
    return 0


def write_metadata_to_bq(metadata, bq_dataset, bq_table):
    """Wrapper function to gather metadata and write to BigQuery
    """
    bq_filename = str(metadata.study_file_id) + '.json'
    push_status = push_metadata_to_bq(metadata, bq_filename, bq_dataset, bq_table)
    return push_status


def check_if_old_output():
    """Exit if old output files found
    """
    output_files = [
        'scp_validation_errors.txt',
        'scp_validation_warnings.txt',
        'bq.json',
    ]

    old_output = False
    for file in output_files:
        if os.path.exists(file):
            print(f'{file} already exists, please delete file and try again')
            old_output = True
    if old_output:
        exit(1)


if __name__ == '__main__':
    args = create_parser().parse_args()
    arguments = vars(args)
    if not args.force:
        check_if_old_output()

    with open(args.convention, 'r') as f:
        convention = json.load(f)
    metadata = CellMetadata(
        file_path=args.input_metadata,
        study_id=args.study_id,
        study_file_id=args.study_file_id,
        study_accession=args.study_accession,
    )
    print('Validating', args.input_metadata)

    validate_input_metadata(metadata, convention, args.bq_json)
    if args.issues_json:
        serialize_issues(metadata)
    report_issues(metadata)
    if args.upload:
        write_metadata_to_bq(metadata, args.bq_dataset, args.bq_table)
    exit_if_errors(metadata)
