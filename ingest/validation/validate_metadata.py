"""Validate input metadata TSV file against metadata convention.

DESCRIPTION
This CLI takes a TSV metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using JSON file for latest Alexandria metadata convention in repo, validate input TSV
$ python3 validate_metadata.py  ../../tests/data/valid_no_array_v2.0.0.tsv

# generate an issues.json file to compare with reference test files
$ python3 validate_metadata.py --issues-json ../../tests/data/valid_no_array_v2.0.0.tsv

# generate a BigQuery upload file to compare with reference test files
$ python3 validate_metadata.py --bq-json ../../tests/data/valid_no_array_v2.0.0.tsv

# use a different metadata convention for validation
$ python3 validate_metadata.py --convention <path to convention json> ../../tests/data/valid_no_array_v2.0.0.tsv

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
import csv
import copy
import itertools
import math

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
    from ..monitor import setup_logger


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
    # make bogus defaults obviously artificial for ease of detection
    parser.add_argument(
        '--study-id',
        help='MongoDB study identifier',
        default='dec0dedfeed1111111111111',
    )
    parser.add_argument(
        '--study-file-id',
        help='MongoDB file identifier',
        default='addedfeed000000000000000',
    )
    parser.add_argument(
        '--study-accession', help='SCP study accession', default='SCPtest'
    )
    parser.add_argument(
        '--bq-dataset', help='BigQuery dataset identifier', default='cell_metadata'
    )
    parser.add_argument(
        '--bq-table', help='BigQuery table identifier', default='alexandria_convention'
    )
    parser.add_argument(
        '--convention',
        help='Metadata convention JSON file',
        default='../../schema/alexandria_convention/alexandria_convention_schema.json',
    )
    parser.add_argument('input_metadata', help='Metadata TSV file')
    return parser


######################## ONTOLOGY RETRIVER  #########################
# TODO: move code in this section to a separate file


# Configure maximum number of seconds to spend & total attempts at external HTTP requests to services, e.g. OLS
MAX_HTTP_REQUEST_TIME = 120
MAX_HTTP_ATTEMPTS = 8


def backoff_handler(details):
    """Handler function to log backoff attempts when querying OLS
    """
    info_logger.debug(
        "Backing off {wait:0.1f} seconds after {tries} tries "
        "calling function {target} with args {args} and kwargs "
        "{kwargs}".format(**details),
        extra={'study_id': None, 'duration': None},
    )


# contains methods for looking up terms in various ontologies,
# as well as caching results of previous queries to speed up performance
class OntologyRetriever:
    cached_ontologies = {}
    cached_terms = {}

    def retrieve_ontology_term_label(self, term, property_name, convention):
        """Retrieve an individual term label from an ontology
        returns JSON payload of ontology, or None if unsuccessful
        Will store any retrieved terms for faster validation of downstream terms

        throws ValueError if the term does not exist or is malformatted
        """

        # initialize cached terms for this property if we haven't seen it yet
        if property_name not in self.cached_terms:
            self.cached_terms[property_name] = {}

        if term not in self.cached_terms[property_name]:
            ontology_urls = get_urls_for_property(convention, property_name)
            self.cached_terms[property_name][
                term
            ] = self.retrieve_ontology_term_label_remote(
                term, property_name, ontology_urls
            )

        return self.cached_terms[property_name][term]

    def retrieve_ontology_term_label_remote(self, term, property_name, ontology_urls):
        """Look up official ontology term from an id, always checking against the remote
        """
        if property_name == 'organ_region':
            return self.retrieve_mouse_brain_term(term, property_name)
        else:
            return self.retrieve_ols_term(ontology_urls, term, property_name)

    # Attach exponential backoff to external HTTP requests
    @backoff.on_exception(
        backoff.expo,
        requests.exceptions.RequestException,
        max_time=MAX_HTTP_REQUEST_TIME,
        max_tries=MAX_HTTP_ATTEMPTS,
        on_backoff=backoff_handler,
        logger="info_logger",
    )
    def retrieve_ols_term(self, ontology_urls, term, property_name):
        """Retrieve an individual term from an ontology
        returns JSON payload of ontology, or None if unsuccessful
        Will store any retrieved ontologies for faster validation of downstream terms
        """
        OLS_BASE_URL = 'https://www.ebi.ac.uk/ols/api/ontologies/'
        # separate ontology shortname from term ID number
        # valid separators are underscore and colon (used by HCA)
        try:
            ontology_shortname, term_id = re.split('[_:]', term)
        except (ValueError, TypeError):
            raise ValueError(
                f'{property_name}: Could not parse provided ontology id, "{term}"'
            )
        # check if we have already retrieved this ontology reference
        if ontology_shortname not in self.cached_ontologies:
            metadata_url = OLS_BASE_URL + ontology_shortname
            self.cached_ontologies[ontology_shortname] = request_json_with_backoff(
                metadata_url
            )
        metadata_ontology = self.cached_ontologies[ontology_shortname]

        # check if the ontology parsed from the term is the same ontology defined in the convention
        # if so, skip the extra call to OLS; otherwise, retrieve the convention-defined ontology for term lookups
        convention_ontology = None
        if metadata_ontology is not None:
            reference_url = metadata_ontology['_links']['self']['href']
            matches = [
                url for url in ontology_urls if url.lower() == reference_url.lower()
            ]
            if matches:
                convention_ontology = metadata_ontology.copy()
            else:
                # store all convention ontologies for lookup later
                for convention_url in ontology_urls:
                    convention_shortname = extract_terminal_pathname(convention_url)
                    convention_ontology = request_json_with_backoff(convention_url)
                    self.cached_ontologies[convention_shortname] = convention_ontology

        if convention_ontology and metadata_ontology:
            base_term_uri = metadata_ontology['config']['baseUris'][0]
            query_iri = encode_term_iri(term_id, base_term_uri)

            term_url = convention_ontology['_links']['terms']['href'] + '/' + query_iri
            # add timeout to prevent request from hanging indefinitely
            response = requests.get(term_url, timeout=60)
            # inserting sleep to minimize 'Connection timed out' error with too many concurrent requests
            time.sleep(0.25)
            if response.status_code == 200:
                return response.json()['label']
            else:
                error_msg = f'{property_name}: No match found in EBI OLS for provided ontology ID: {term}'
                raise ValueError(error_msg)
        elif not metadata_ontology:
            error_msg = f'No result from EBI OLS for provided ontology shortname \"{ontology_shortname}\"'
            print(error_msg)
            info_logger.info(error_msg, extra={'study_id': None, 'duration': None})
            raise ValueError(
                f'{property_name}: No match found in EBI OLS for provided ontology ID: {term}'
            )
        else:
            error_msg = (
                f'encountered issue retrieving {ontology_urls} or {ontology_shortname}'
            )
            print(error_msg)
            info_logger.info(error_msg, extra={'study_id': None, 'duration': None})
            raise RuntimeError(error_msg)

    def retrieve_mouse_brain_term(self, term, property_name):
        mouse_brain_atlas = self.fetch_allen_mouse_brain_atlas()
        MBA_id = parse_organ_region_ontology_id(term)
        if MBA_id not in mouse_brain_atlas:
            raise ValueError(
                f'{property_name}: No match found in Allen Mouse Brain Atlas for provided ontology ID: {term}'
            )
        return mouse_brain_atlas[MBA_id]

    def fetch_allen_mouse_brain_atlas(self):
        """Get the mouse brain atlas file, either remote or cached, and parse it into an id -> name dictionary
        """
        if 'allen_mouse_brain_atlas' not in self.cached_ontologies:
            self.cached_ontologies[
                'allen_mouse_brain_atlas'
            ] = fetch_allen_mouse_brain_atlas_remote()
        return self.cached_ontologies['allen_mouse_brain_atlas']


def parse_organ_region_ontology_id(term):
    """Extract term id from valid identifiers or raise a ValueError
    """
    try:
        ontology_shortname, term_id = re.split('[_:]', term)
        if ontology_shortname == 'MBA':
            # remove leading zeroes (e.g. '000786' => '786') since the dictionary file doesn't have them
            return term_id.lstrip('0')
        else:
            error_msg = f'organ_region: Invalid ontology code, "{ontology_shortname}"'
            raise ValueError(error_msg)
    except (TypeError, ValueError):
        # when term value is empty string -> TypeError, convert this to a value error
        raise ValueError(
            f'organ_region: Could not parse provided ontology id, "{term}"'
        )


def fetch_allen_mouse_brain_atlas_remote():
    """Get the mouse brain atlas file from the remote source, and parse it into an id -> name dictionary
    """
    MBA_url = 'https://raw.githubusercontent.com/broadinstitute/scp-ingest-pipeline/58f9308e425c667a34219a3dcadf7209fe09f788/schema/organ_region/mouse_brain_atlas/MouseBrainAtlas.csv'
    region_dict = {}
    try:
        region_ontology = requests.get(MBA_url).text
        csv_data = csv.DictReader(region_ontology.splitlines())
        for row in csv_data:
            region_dict[row['id']] = row['name']
        return region_dict
    # if we can't get the file in GitHub for validation record error
    except:  # noqa E722
        error_msg = 'Unable to read GitHub-hosted organ_region ontology file'
        raise RuntimeError(error_msg)


# Attach exponential backoff to external HTTP requests
@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=MAX_HTTP_REQUEST_TIME,
    max_tries=MAX_HTTP_ATTEMPTS,
    on_backoff=backoff_handler,
    logger="info_logger",
)
def request_json_with_backoff(ontology_url):
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


def get_urls_for_property(convention, property_name):
    return convention['properties'][property_name]['ontology'].split(',')


######################## END ONTOLOGY RETRIVER DEFINITION #########################

# create an OntologyRetriever instance to handle fetching and caching ontology terms
retriever = OntologyRetriever()


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


def is_required_metadata(convention, metadatum):
    """Check in metadata convention if metadata is required
    """
    try:
        if metadatum in convention['required']:
            return True
        else:
            return False
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


def insert_array_ontology_label_row_data(
    property_name, row, metadata, required, convention, ontology_label
):
    cell_id = row['CellID']
    # if there are differing numbers of id/label values, and not empty, log error and don't try to fix
    if len(row[property_name]) != len(row[ontology_label]) and row[ontology_label]:
        error_msg = f'{property_name}: mismatched # of {property_name} and {ontology_label} values'
        metadata.store_validation_issue('error', 'ontology', error_msg, [cell_id])
        return row
    if not row[ontology_label]:
        array_label_for_bq = []
        for id in row[property_name]:
            label_lookup = ''
            try:
                label_lookup = retriever.retrieve_ontology_term_label(
                    id, property_name, convention
                )
                reference_ontology = (
                    "EBI OLS lookup"
                    if property_name != "organ_region"
                    else "Mouse Brain Atlas ontology"
                )
                error_msg = (
                    f'{property_name}: missing ontology label '
                    f'"{id}" - using "{label_lookup}" per {reference_ontology}'
                )
                metadata.store_validation_issue(
                    'warn', 'ontology', error_msg, [cell_id]
                )
            except BaseException as e:
                print(e)
                error_msg = (
                    f'{property_name}: unable to lookup "{id}" - '
                    f'cannot propagate array of labels for {property_name};'
                    f'may cause inaccurate error reporting for {property_name}'
                )
                metadata.store_validation_issue(
                    'error', 'ontology', error_msg, [cell_id]
                )
            array_label_for_bq.append(label_lookup)
        row[ontology_label] = array_label_for_bq

    for id, label in zip(row[property_name], row[ontology_label]):
        metadata.ontology[property_name][(id, label)].append(cell_id)
    return row


def insert_ontology_label_row_data(
    property_name, row, metadata, required, convention, ontology_label
):
    id = row[property_name]
    cell_id = row['CellID']

    if not row[ontology_label]:
        # for optional columns, try to fill it in
        try:
            label = retriever.retrieve_ontology_term_label(
                id, property_name, convention
            )
            row[ontology_label] = label
            reference_ontology = (
                "EBI OLS lookup"
                if property_name != "organ_region"
                else "Mouse Brain Atlas ontology"
            )
            error_msg = (
                f'{property_name}: missing ontology label '
                f'"{id}" - using "{label}" per {reference_ontology}'
            )
            metadata.store_validation_issue('warn', 'ontology', error_msg, [cell_id])
        except BaseException as e:
            print(e)
            error_msg = f'Optional column {ontology_label} empty and could not be resolved from {property_name} column value {row[property_name]}'
            metadata.store_validation_issue('warn', 'ontology', error_msg, [cell_id])

    metadata.ontology[property_name][(row[property_name], row[ontology_label])].append(
        cell_id
    )
    return row


def collect_cell_for_ontology(
    property_name, row, metadata, convention, array, required
):
    """Collect ontology info for a single metadatum into CellMetadata.ontology dictionary
    if ontology_label not provided for optional metadata, label is looked up
    and added to BigQuery json to populate metadata backend
    Missing ontology_label for required metadata is not inserted, should fail validation
    """

    if property_name.endswith('__unit'):
        ontology_label = property_name + '_label'
    else:
        ontology_label = property_name + '__ontology_label'
    updated_row = copy.deepcopy(row)
    cell_id = updated_row['CellID']

    # for case where they've omitted a column altogether or left it blank, add a blank entry
    if ontology_label not in updated_row or value_is_nan(updated_row[ontology_label]):
        if array:
            updated_row[ontology_label] = []
        else:
            updated_row[ontology_label] = ''

    # !!!property_name checking is probably not needed
    if property_name not in updated_row or value_is_nan(updated_row[property_name]):
        if array:
            updated_row[property_name] = []
        else:
            updated_row[property_name] = ''

    if (not updated_row[ontology_label] or not updated_row[property_name]) and required:
        # Catch cases where ontology_label or property_name column completely empty
        missing_column_message = ''
        if not updated_row[property_name]:
            missing_column_message = (
                f'{property_name}: required column "{property_name}" missing data'
            )
        if not updated_row[ontology_label]:
            if not updated_row[property_name]:
                missing_column_message += ' and '
            missing_column_message += (
                f'{property_name}: required column "{ontology_label}" missing data'
            )
        # for required columns, just log the error and continue
        metadata.store_validation_issue(
            'error', 'ontology', missing_column_message, [cell_id]
        )
        # metadata.ontology[property_name][(updated_row[property_name], None)].append(cell_id)
    else:
        if array:
            updated_row = insert_array_ontology_label_row_data(
                property_name,
                updated_row,
                metadata,
                required,
                convention,
                ontology_label,
            )
        else:
            updated_row = insert_ontology_label_row_data(
                property_name,
                updated_row,
                metadata,
                required,
                convention,
                ontology_label,
            )

    return updated_row


def collect_ontology_data(row_data, metadata, convention):
    """Collect unique ontology IDs for ontology validation
    """
    row_ref = copy.deepcopy(row_data)
    for entry in row_ref.keys():
        if is_ontology_metadata(convention, entry):
            updated_row = collect_cell_for_ontology(
                entry,
                row_data,
                metadata,
                convention,
                array=is_array_metadata(convention, entry),
                required=is_required_metadata(convention, entry),
            )
            # update BigQuery data to include missing ontology labels
            # for optional ontology metadata
            row_data.update(updated_row)
    return row_data


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
    metadata_annots = dict(itertools.zip_longest(metadata_names, type_annots))
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


def value_is_nan(value):
    """Check if value is nan
    nan is a special dataframe value to indicate missing data
    """
    try:
        return math.isnan(value)
    except TypeError:
        return False


def cast_integer_type(value):
    """Cast metadata value as integer
    """
    if value_is_nan(value):
        # nan indicates missing data, has no valid integer value for casting
        return value
    else:
        return int(value)


def cast_float_type(value):
    """Cast metadata value as float
    """
    return float(value)


def cast_string_type(value):
    """Cast string type per convention where Pandas autodetected a number
    """
    if value_is_nan(value):
        # nan indicates missing data; by type, nan is a numpy float
        # so a separate type check is needed for proper handling
        return value
    elif isinstance(value, numbers.Number):
        return str(value)
    else:
        return value


def regularize_ontology_id(value):
    """Regularize ontology_ids for storage with underscore format
    """
    try:
        return value.replace(":", "_")
    except AttributeError:
        # when expected value is not actually an ontology ID
        # return the bad value for JSON schema validation
        return value


def cast_metadata_type(metadatum, value, id_for_error_detail, convention, metadata):
    """For metadatum, lookup expected type by metadata convention
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
                try:
                    if 'ontology' in convention['properties'][metadatum]:
                        element = regularize_ontology_id(element)
                except KeyError:
                    # study-specific metadata (ie. non-metadata-convention metadata)
                    # should no longer trigger this exception but should be allowed to pass
                    pass
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
            try:
                if 'ontology' in convention['properties'][metadatum]:
                    value = regularize_ontology_id(value)
            except KeyError:
                # study-specific metadata (ie. non-metadata convention metadata)
                # should no longer trigger this exception but should be allowed to pass
                pass
            cast_metadata[metadatum] = [
                metadata_types.get(lookup_metadata_type(convention, metadatum))(value)
            ]
    else:
        try:
            if 'ontology' in convention['properties'][metadatum]:
                value = regularize_ontology_id(value)
        except KeyError:
            # study-specific metadata (ie. non-metadata convention metadata)
            # should no longer trigger this exception but should be allowed to pass
            pass
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
    check metadata type is of type assigned in convention
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
    row_info = dict(itertools.zip_longest(metadata_names, line))
    processed_row = {}
    for k, v in row_info.items():
        # for metadata not in convention, no need to process
        if k not in convention['properties'].keys():
            continue
        # for optional metadata, do not pass empty cells (nan)
        if k not in convention['required'] and value_is_nan(v):
            continue
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
            row = collect_ontology_data(row, metadata, convention)
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
                # extract convention version from URI in convention JSON file
                row['metadata_convention_version'] = convention['$id'].split("/")[-2]
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
    """Print issue to console with coloring and
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


def validate_collected_ontology_data(metadata, convention):
    """Evaluate collected ontology_id, ontology_label info in
    CellMetadata.ontology dictionary by querying ontology source
    validity of ontology_id and cross-check that input ontology_label
    matches ontology_label from ontology lookup
    """
    for property_name in metadata.ontology.keys():
        # pull in file for validation of non-EBI OLS metadata for organ_region
        if property_name == 'organ_region':
            try:
                region_ontology = (  # noqa: F841
                    retriever.fetch_allen_mouse_brain_atlas()
                )
            except BaseException as e:
                print(e)
                # immediately return as validation should not continue if ontology file unavailable
                return
        # split on comma in case this property from the convention supports multiple ontologies
        ontology_urls = convention['properties'][property_name]['ontology'].split(',')

        for ontology_info in metadata.ontology[property_name].keys():
            ontology_id, ontology_label = ontology_info
            try:
                matched_label_for_id = retriever.retrieve_ontology_term_label(
                    ontology_id, property_name, convention
                )
                if matched_label_for_id != ontology_label:
                    ontology_source_name = 'EBI OLS'
                    if property_name == 'organ_region':
                        ontology_source_name = 'Allen Mouse Brain Atlas'
                    if is_required_metadata(convention, property_name) and (
                        not ontology_label or value_is_nan(ontology_label)
                    ):
                        error_msg = (
                            f'{property_name}: ontology label required for \"{ontology_id}\" '
                            f'to cross-check for data entry error"'
                        )
                        metadata.store_validation_issue(
                            'error',
                            'ontology',
                            error_msg,
                            metadata.ontology[property_name][
                                (ontology_id, ontology_label)
                            ],
                        )
                    else:
                        error_msg = (
                            f'{property_name}: input ontology_label \"{ontology_label}\" '
                            f'does not match {ontology_source_name} lookup \"{matched_label_for_id}\" for ontology id \"{ontology_id}\"'
                        )
                        metadata.store_validation_issue(
                            'error',
                            'ontology',
                            error_msg,
                            metadata.ontology[property_name][
                                (ontology_id, ontology_label)
                            ],
                        )
            except ValueError as valueError:
                metadata.store_validation_issue(
                    'error',
                    'ontology',
                    valueError.args[0],
                    metadata.ontology[property_name][(ontology_id, ontology_label)],
                )
            except requests.exceptions.RequestException as err:
                error_msg = f'External service outage connecting to {ontology_urls} when querying {ontology_id}:{ontology_label}: {err}'
                error_logger.error(error_msg)
                metadata.store_validation_issue('error', 'ontology', error_msg),
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
            # existence of unit metadata is enforced by jsonschema
            # any empty cells in a unit column are okay to be empty
            if metadata.file[name].nunique(dropna=True).values[0] != 1:
                error_msg = (
                    f'{name}: values for each unit metadata required to be uniform'
                )
                metadata.store_validation_issue('error', 'convention', error_msg)


def serialize_bq(bq_dict, filename='bq.json'):
    """Write metadata collected for validation to json file
    BigQuery requires newline delimited json objects
    """
    bq_dict_copy = bq_dict.copy()
    if 'organism_age' in bq_dict and 'organism_age__unit_label' in bq_dict:
        bq_dict_copy['organism_age__seconds'] = calculate_organism_age_in_seconds(
            bq_dict['organism_age'], bq_dict['organism_age__unit_label']
        )

    data = json.dumps(bq_dict_copy)
    with open(filename, 'a') as jsonfile:
        jsonfile.write(data + '\n')


def calculate_organism_age_in_seconds(organism_age, organism_age_unit_label):
    multipliers = {
        'microsecond': 0.000001,
        'millisecond': 0.001,
        'second': 1,
        'minute': 60,
        'hour': 3600,
        'day': 86400,
        'week': 604800,
        'month': 2626560,  # (day * 30.4 to fuzzy-account for different months)
        'year': 31557600,  # (day * 365.25 to fuzzy-account for leap-years)
    }
    multiplier = multipliers[organism_age_unit_label]
    return organism_age * multiplier


def serialize_issues(metadata):
    """Write collected issues to json file
    """
    with open('issues.json', 'w') as jsonfile:
        json.dump(metadata.issues, jsonfile, indent=2)
        jsonfile.write("\n")  # Add newline cause Py JSON does not


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
    """Upload local NDJSON to BigQuery
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
