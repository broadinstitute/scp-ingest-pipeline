"""Validate input metadata TSV file against metadata convention.
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
import pandas as pd
import gzip
import glob

import colorama
from colorama import Fore
import jsonschema
from google.cloud import bigquery

sys.path.append("..")
try:
    # Used when importing internally and in tests
    from monitor import setup_logger
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..monitor import setup_logger


dev_logger = setup_logger(__name__, "log.txt", format="support_configs")

user_logger = setup_logger(
    __name__ + ".user_logger", "user_log.txt", level=logging.ERROR
)

# Ensures normal color for print() output, unless explicitly changed
colorama.init(autoreset=True)

# Configure maximum number of seconds to spend & total attempts at external HTTP requests to services, e.g. OLS
MAX_HTTP_REQUEST_TIME = 120
MAX_HTTP_ATTEMPTS = 8


######################## ONTOLOGY RETRIVER  #########################
# TODO: move code in this section to a separate file


# Configure maximum number of seconds to spend & total attempts at external HTTP requests to services, e.g. OLS
MAX_HTTP_REQUEST_TIME = 120
MAX_HTTP_ATTEMPTS = 8


def backoff_handler(details):
    """Handler function to log backoff attempts when querying OLS"""
    dev_logger.debug(
        "Backing off {wait:0.1f} seconds after {tries} tries "
        "calling function {target} with args {args} and kwargs "
        "{kwargs}".format(**details)
    )

# handles reading minified ontologies and performing term/synonym lookups
class MinifiedOntologyReader():
    parsed_ontologies = {}

    def __init__(self):
        ontology_dir = f"{os.path.dirname(os.path.realpath(__file__))}/ontologies"
        for ontology_file in glob.glob(f"{ontology_dir}/*.min.tsv.gz"):
            ontology_name = ontology_file.split('/')[-1].replace(".min.tsv.gz", "")
            self.populate_ontology(ontology_name, ontology_file)

    def ontology_names(self):
        return list(self.parsed_ontologies.keys())

    def populate_ontology(self, ontology_name, ontology_file):
        """Parses ontology file by name and populates entries into parsed_ontologies for lookup
        :param ontology_name: name of ontology
        :param ontology_file: relative path to ontology file
        :return: parsed ontology dictionary
        """
        dev_logger.debug(f"populating minified ontology {ontology_name} from {ontology_file}")
        with gzip.open(ontology_file, 'rt') as file_gz:
            ontology = {}
            for line in file_gz.readlines():
                try:
                    ontology_id, label, raw_syn = line.split("\t")
                    entry = {"label": label, "synonyms": [syn.replace("\n", '') for syn in raw_syn.split("||")]}
                    ontology[ontology_id] = entry
                except (ValueError, TypeError) as e:
                    dev_logger.error(f"could not process {line} from {ontology_name}: {e}")
            self.parsed_ontologies[ontology_name] = ontology

    def find_ontology_entry(self, ontology_name, identifier, property_name):
        """Find an entry in a parsed ontology by identfier
        :param ontology_name: name of ontology
        :param identifier: ontology ID, e.g. MONDO_0005887
        :param property_name: name of metadata property, e.g. species
        :return: dict
        """
        entry = self.parsed_ontologies.get(ontology_name, {}).get(identifier, {})
        if entry:
            return entry
        else:

            msg = f"{property_name}: No match found in EBI OLS for provided ontology ID: {identifier}"
            raise ValueError(msg)


# contains methods for looking up terms in various ontologies,
# as well as caching results of previous queries to speed up performance
class OntologyRetriever:
    cached_ontologies = {}
    cached_terms = {}

    def retrieve_ontology_term_label_and_synonyms(
        self, term, property_name, convention, attribute_type: str
    ):
        """Retrieve an individual term label and any synonymns from an ontology
        returns JSON payload of ontology, or None if unsuccessful
        Will store any retrieved labels/synonyms for faster validation of downstream terms

        throws ValueError if the term does not exist or is malformatted
                :param term: ontology term
                :param property_name: term attribute Ex. biosample_id, CellID, etc.
                :param convention: metadata convention being checked against
                :param attribute_type: attribute type for term (string, array, boolean)
        """

        # initialize cached terms for this property if we haven't seen it yet
        if property_name not in self.cached_terms:
            self.cached_terms[property_name] = {}

        if term not in self.cached_terms[property_name]:
            ontology_urls = get_urls_for_property(convention, property_name)
            self.cached_terms[property_name][term] = (
                self.retrieve_ontology_term_label_remote(
                    term, property_name, ontology_urls, attribute_type
                )
            )

        return self.cached_terms[property_name][term]

    def retrieve_ontology_term_label_remote(
        self, term, property_name, ontology_urls, attribute_type
    ):
        """Look up official ontology term from an id, always checking against the remote"""
        # organ_region currently uses single a non-OLS ontology, the Allen Institute's Mouse Brain Atlas (MBA)
        # MBA was originally downloaded April 2020 from Allen Brain Atlas data portal URL
        # http://api.brain-map.org/api/v2/data/query.csv?criteria=model::Structure,rma::criteria,[ontology_id$eq1],rma::options[order$eq%27structures.graph_order%27][num_rows$eqall]
        # Mouse Brain Atlas IDs are numeric, to make the IDs easily attributable to their ontology
        # SCP study owners are expected to prepend "MBA" and pad IDs shorter than 9 digits with leading zeros
        # the expected formatting brings the MBA IDs into a format similar to other ontologies
        # and avoids potential name collisions
        if property_name == "organ_region":
            return self.retrieve_mouse_brain_term(term, property_name)
        else:
            # leave debug statement for QA purposes later
            dev_logger.debug(
                f"Using fallback EBI OLS call with {ontology_urls}, {term}, {property_name}"
            )
            return self.retrieve_ols_term(
                ontology_urls, term, property_name, attribute_type
            )

    # Attach exponential backoff to external HTTP requests
    @backoff.on_exception(
        backoff.expo,
        requests.exceptions.RequestException,
        max_time=MAX_HTTP_REQUEST_TIME,
        max_tries=MAX_HTTP_ATTEMPTS,
        on_backoff=backoff_handler,
        logger="dev_logger",
    )
    def retrieve_ols_term(self, ontology_urls, term, property_name, attribute_type):
        """Retrieve an individual term from an ontology
        returns JSON payload of ontology, or None if unsuccessful
        Will store any retrieved ontologies for faster validation of downstream terms
        """
        OLS_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/"
        # separate ontology shortname from term ID number
        # valid separators are underscore and colon (used by HCA)
        try:
            ontology_shortname, term_id = re.split("[_:]", term)
        except (ValueError, TypeError):
            msg = f'{property_name}: Could not parse provided ontology id, "{term}".'
            raise ValueError(msg)
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
            reference_url = metadata_ontology["_links"]["self"]["href"]
            matches = [
                url for url in ontology_urls if url.lower() == reference_url.lower()
            ]
            if matches:
                convention_ontology = metadata_ontology.copy()
            else:
                # Store all convention ontologies for lookup later.
                # Prioritize first ontology, if multiple as with MONDO, PATO in disease
                for convention_url in list(reversed(ontology_urls)):
                    convention_shortname = extract_terminal_pathname(convention_url)
                    convention_ontology = request_json_with_backoff(convention_url)
                    self.cached_ontologies[convention_shortname] = convention_ontology

        if convention_ontology and metadata_ontology:
            base_term_uri = get_ontology_file_location(metadata_ontology)
            query_iri = encode_term_iri(term, base_term_uri)

            term_url = convention_ontology["_links"]["terms"]["href"] + "/" + query_iri

            # add timeout to prevent request from hanging indefinitely
            response = requests.get(term_url, timeout=60)
            # inserting sleep to minimize 'Connection timed out' error with too many concurrent requests
            time.sleep(0.25)
            if response.status_code == 200:
                # return canonical label, along with synonyms and any 'related synonyms' to match later
                term_json = response.json()
                labels = {"label": term_json["label"], "synonyms": []}
                synonyms = []
                if term_json['synonyms']:
                    synonyms += term_json['synonyms']
                # safe lookup of nested dictionary
                related_synonyms = term_json.get('annotation', {}).get(
                    'has_related_synonym'
                )
                if related_synonyms:
                    synonyms += related_synonyms
                # safe lookup of nested dictionary
                exact_synonyms = term_json.get('annotation', {}).get(
                    'has_exact_synonym'
                )
                if exact_synonyms:
                    synonyms += exact_synonyms
                # uniquify list via set and return
                labels["synonyms"] = list(set(synonyms))
                return labels
            else:
                msg = f"{property_name}: No match found in EBI OLS for provided ontology ID: {term}"
                raise ValueError(msg)
        elif not metadata_ontology:
            msg = f'No result from EBI OLS for provided ontology shortname "{ontology_shortname}".'
            print(msg)
            user_logger.error(msg)
            raise ValueError(
                f"{property_name}: No match found in EBI OLS for provided ontology ID: {term}"
            )
        else:
            msg = (
                f"Encountered issue retrieving {ontology_urls} or {ontology_shortname}."
            )
            print(msg)
            user_logger.info(msg)
            raise RuntimeError(msg)

    def retrieve_mouse_brain_term(self, term, property_name):
        """Determine whether ID is in mouse brain atlas (MBA) file"""
        # MBA ID is also the leaf entity of structure_id_path in the MBA file
        # Entries with short structure_id_path seem to be synonymous with
        # Uberon terms, suggesting MBA terms could be mapped as extensions of the
        # Uberon ontology in future
        mouse_brain_atlas = self.fetch_allen_mouse_brain_atlas()
        MBA_id = parse_organ_region_ontology_id(term)
        if MBA_id not in mouse_brain_atlas:
            raise ValueError(
                f"{property_name}: No match found in Allen Mouse Brain Atlas for provided ontology ID: {term}"
            )
        return {"label": mouse_brain_atlas[MBA_id]}

    def fetch_allen_mouse_brain_atlas(self):
        """Get the mouse brain atlas file, either remote or cached, and parse it into an id -> name dictionary"""
        if "allen_mouse_brain_atlas" not in self.cached_ontologies:
            self.cached_ontologies["allen_mouse_brain_atlas"] = (
                fetch_allen_mouse_brain_atlas_remote()
            )
        return self.cached_ontologies["allen_mouse_brain_atlas"]


def parse_organ_region_ontology_id(term):
    """Extract term id from valid identifiers or raise a ValueError"""
    try:
        ontology_shortname, term_id = re.split("[_:]", term)
        if ontology_shortname == "MBA":
            # remove leading zeroes (e.g. '000786' => '786') since the dictionary file doesn't have them
            return term_id.lstrip("0")
        else:
            msg = f'organ_region: Invalid ontology code, "{ontology_shortname}".'
            raise ValueError(msg)
    except (TypeError, ValueError):
        # when term value is empty string -> TypeError, convert this to a value error
        raise ValueError(
            f'organ_region: Could not parse provided ontology ID, "{term}".'
        )


def fetch_allen_mouse_brain_atlas_remote():
    """Get the mouse brain atlas file from the remote source, and parse it into an id -> name dictionary"""
    MBA_url = "https://raw.githubusercontent.com/broadinstitute/scp-ingest-pipeline/58f9308e425c667a34219a3dcadf7209fe09f788/schema/organ_region/mouse_brain_atlas/MouseBrainAtlas.csv"
    region_dict = {}
    try:
        region_ontology = requests.get(MBA_url).text
        csv_data = csv.DictReader(region_ontology.splitlines())
        for row in csv_data:
            region_dict[row["id"]] = row["name"]
        return region_dict
    # if we can't get the file in GitHub for validation record error
    except:  # noqa E722
        msg = "Unable to read GitHub-hosted organ_region ontology file."
        dev_logger.exception(msg)
        raise RuntimeError(msg)


# Attach exponential backoff to external HTTP requests
@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=MAX_HTTP_REQUEST_TIME,
    max_tries=MAX_HTTP_ATTEMPTS,
    on_backoff=backoff_handler,
    logger="dev_logger",
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
    """Extract the last path segment from a URL"""
    return list(filter(None, url.split("/"))).pop()


def get_urls_for_property(convention, property_name):
    return convention["properties"][property_name]["ontology"].split(",")


def get_ontology_file_location(ontology):
    """Find and format fileLocation URL for running queries against an ontology"""
    location = ontology["config"]["fileLocation"]
    # issue with some ontologies (like GO) having non-standard fileLocation URLs
    if 'extensions' in location:
        location = location.replace(f"{ontology['ontologyId']}/extensions/", "")
    if location == "http://purl.obolibrary.org/obo/NCBITAXON_":
        location = "http://purl.obolibrary.org/obo/NCBITaxon_"
    return '/'.join(location.split('/')[:-1]) + '/'


######################## END ONTOLOGY RETRIVER DEFINITION #########################

# create an OntologyRetriever instance to handle fetching and caching ontology terms
retriever = OntologyRetriever()
minified_reader = MinifiedOntologyReader()

def validate_schema(json, metadata):
    """Check validity of metadata convention as JSON schema.
    if valid, return jsonschema validator object else return None
    """

    try:
        jsonschema.Draft7Validator.check_schema(json)
        valid_schema = jsonschema.Draft7Validator(json)
        return valid_schema
    except jsonschema.SchemaError:
        msg = "Invalid metadata convention file, cannot validate metadata."
        metadata.store_validation_issue(
            "error", msg, "runtime:invalid-convention", issue_type="runtime"
        )
        return None


def is_array_metadata(convention, metadatum):
    """Check if metadata is array type from metadata convention"""
    try:
        type_lookup = convention["properties"][metadatum]["type"]
        if type_lookup == "array":
            return True
        else:
            return False
    except KeyError:
        return False


def is_ontology_metadata(convention, metadatum):
    """Check if metadata is ontology from metadata convention"""
    try:
        return bool(convention["properties"][metadatum]["ontology"])
    except KeyError:
        return False


def is_required_metadata(convention, metadatum):
    """Check in metadata convention if metadata is required"""
    try:
        if metadatum in convention["required"]:
            return True
        else:
            return False
    except KeyError:
        return False


def lookup_metadata_type(convention, metadatum):
    """Look up metadata type from metadata convention"""
    try:
        if is_array_metadata(convention, metadatum):
            type_lookup = convention["properties"][metadatum]["items"]["type"]
        else:
            type_lookup = convention["properties"][metadatum]["type"]
        return type_lookup
    except KeyError:
        return None


def list_duplicates(cells):
    """Find duplicates in list for detailed reporting"""
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
    """Check all CellID are unique."""
    valid = False
    if len(metadata.cells) == len(set(metadata.cells)):
        valid = True
    else:
        dups = list_duplicates(metadata.cells)
        msg = "Duplicate CellID(s) in metadata file."
        metadata.store_validation_issue(
            "error", msg, "content:duplicate:cells-within-file", associated_info=dups
        )
    return valid

def retrieve_label_and_synonyms(
    ontology_id, property_name, convention, property_type
):
    """Wrapper method to retrieve label and synonyms depending on whether ontology is local or remote
    :param ontology_id: ontology ID, e.g. MONDO_0005887
    :param property_name: name of metadata property, e.g. species
    :param convention: metadata convention being checked against
    :param property_type: attribute type for term (string, array, boolean)
    """
    ontology_name = re.split("[_:]", ontology_id)[0].lower()
    if ontology_is_local(ontology_name):
        return minified_reader.find_ontology_entry(ontology_name, ontology_id, property_name)
    else:
        return retriever.retrieve_ontology_term_label_and_synonyms(
            ontology_id, property_name, convention, property_type
        )

def insert_array_ontology_label_row_data(
    property_name, row, metadata, required, convention, ontology_label
):
    cell_id = row["CellID"]
    # if there are differing numbers of id/label values, and not empty, log error and don't try to fix
    if len(row[property_name]) != len(row[ontology_label]) and row[ontology_label]:
        msg = f"{property_name}: mismatched # of {property_name} and {ontology_label} values."
        metadata.store_validation_issue(
            "error", msg, "ontology:array-length-mismatch", associated_info=[cell_id]
        )
        return row

    if not row[ontology_label]:
        array_label_for_bq = []
        # track original labels, including blanks, in the ordered ontology structure
        metadata.ordered_ontology[property_name].extend(row[property_name])
        metadata.ordered_labels[property_name].extend('')
        for id in row[property_name]:
            label_lookup = ""
            try:

                label_and_synonyms = retrieve_label_and_synonyms(id, property_name, convention, "array")
                label_lookup = label_and_synonyms.get('label')
                reference_ontology = (
                    "EBI OLS lookup"
                    if property_name != "organ_region"
                    else "Mouse Brain Atlas ontology"
                )
                msg = (
                    f"{property_name}: missing ontology label "
                    f'"{id}" - using "{label_lookup}" per {reference_ontology}'
                )
                metadata.store_validation_issue(
                    "warn",
                    msg,
                    "ontology:missing-label-lookup",
                    associated_info=[cell_id],
                )
            except BaseException as e:
                print(e)
                msg = (
                    f'{property_name}: unable to lookup "{id}" - '
                    f"cannot propagate array of labels for {property_name};"
                    f"may cause inaccurate error reporting for {property_name}."
                )
                metadata.store_validation_issue(
                    "error",
                    msg,
                    "ontology:label-lookup-error",
                    associated_info=[cell_id],
                )
            array_label_for_bq.append(label_lookup)
        row[ontology_label] = array_label_for_bq
    else:
        metadata.ordered_ontology[property_name].extend(row[property_name])
        metadata.ordered_labels[property_name].extend(row[ontology_label])
    for id, label in zip(row[property_name], row[ontology_label]):
        metadata.ontology[property_name][(id, label)].append(cell_id)
    return row


def insert_ontology_label_row_data(
    property_name, row, metadata, required, convention, ontology_label
):
    id = row[property_name]
    cell_id = row["CellID"]

    if not row[ontology_label]:
        # track original labels, including blanks, in the ordered ontology structure
        metadata.ordered_ontology[property_name].append(id)
        metadata.ordered_labels[property_name].append('')
        # for optional columns, try to fill it in
        property_type = convention["properties"][property_name]["type"]
        try:
            label_and_synonyms = retrieve_label_and_synonyms(id, property_name, convention, property_type)
            label = label_and_synonyms.get('label')
            row[ontology_label] = label
            reference_ontology = (
                "EBI OLS lookup"
                if property_name != "organ_region"
                else "Mouse Brain Atlas ontology"
            )
            msg = (
                f"{property_name}: missing ontology label "
                f'"{id}" - using "{label}" per {reference_ontology}'
            )
            metadata.store_validation_issue(
                "warn", msg, "ontology:missing-label-lookup", associated_info=[cell_id]
            )
        except BaseException as e:
            print(e)
            msg = f"Optional column {ontology_label} empty and could not be resolved from {property_name} column value {row[property_name]}."
            metadata.store_validation_issue(
                "warn", msg, "ontology:label-lookup-error", associated_info=[cell_id]
            )
    else:
        metadata.ordered_ontology[property_name].append(id)
        metadata.ordered_labels[property_name].append(row[ontology_label])

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

    def omitted_values(ontology_label, cell_name):
        return (
            # check for missing values
            ontology_label not in updated_row
            # check for empty cells
            or value_is_nan(cell_name)
            # check for empty cells
            or is_empty_string(cell_name)
        )

    if property_name.endswith("__unit"):
        ontology_label = property_name + "_label"
    else:
        ontology_label = property_name + "__ontology_label"
    updated_row = copy.deepcopy(row)
    cell_id = updated_row["CellID"]

    # for case where they've omitted a column altogether or left it blank, add a blank entry
    if omitted_values(ontology_label, updated_row.get(ontology_label)):
        if array:
            updated_row[ontology_label] = []
        else:
            updated_row[ontology_label] = ""

    if omitted_values(property_name, updated_row[property_name]):
        if array:
            updated_row[property_name] = []
        else:
            updated_row[property_name] = ""

    if (not updated_row[ontology_label] or not updated_row[property_name]) and required:
        # Catch cases where ontology_label or property_name column completely empty
        missing_column_message = ""
        if not updated_row[property_name]:
            missing_column_message = (
                f'{property_name}: required column "{property_name}" missing data'
            )
        if not updated_row[ontology_label]:
            if not updated_row[property_name]:
                missing_column_message += " and "
            missing_column_message += (
                f'{property_name}: required column "{ontology_label}" missing data'
            )
        # for required columns, just log the error and continue
        metadata.store_validation_issue(
            "error",
            missing_column_message,
            "content:missing-required-values",
            associated_info=[cell_id],
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
    """Collect unique ontology IDs for ontology validation"""
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
    """Check if metadata type annotation is consistent with metadata convention type"""
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    type_annots = metadata.file.columns.get_level_values(1).tolist()
    # if input was TSV metadata file, SCP format requires 'NAME' for the first
    # column which is expected to be CellID, primarily based on loom convention
    # would do conditional (eg. if metadata_names[0].upper() == 'NAME':)
    # but subsequent code expects 'CellID' even if header check fails and first
    # word in header is not some form of 'NAME' and script breaks
    metadata_names[0] = "CellID"
    type_annots[0] = "group"
    metadata_annots = dict(itertools.zip_longest(metadata_names, type_annots))
    annot_equivalents = {
        "numeric": ["number", "integer"],
        "group": ["boolean", "string"],
    }
    for metadatum, annot in metadata_annots.items():
        convention_type = lookup_metadata_type(convention, metadatum)
        try:
            if convention_type and convention_type not in annot_equivalents.get(annot):
                for k, v in annot_equivalents.items():
                    if convention_type in v:
                        expected = k
                msg = (
                    f'{metadatum}: "{annot}" annotation in metadata file conflicts with metadata convention. '
                    f'Convention expects "{expected}" values.'
                )
                metadata.store_validation_issue(
                    "error", msg, "content:type:value-type-mismatch"
                )
        except TypeError:
            for k, v in annot_equivalents.items():
                if convention_type in v:
                    expected = k
            if "." in annot:
                # duplicated metadata header name identified in validate_unique_header
                # detection of invalid type annotation is side effect of Pandas
                pass
            elif "Unnamed" in annot:
                # missing type annotation also detected in validate_against_header_count
                # invalid type annotation is side effect of Pandas
                msg = (
                    f"{metadatum}: missing TYPE annotation in metadata file. "
                    f'Convention expects "{expected}" annotation.'
                )
                metadata.store_validation_issue("error", msg, "format:cap:type")
            else:
                msg = (
                    f'{metadatum}: invalid "{annot}" annotation in metadata file. '
                    f'Convention expects "{expected}" annotation.'
                )
                metadata.store_validation_issue("error", msg, "format:cap:type")


def cast_boolean_type(value):
    """Cast metadata value as boolean, if castable"""
    if isinstance(value, bool):
        return value
    elif str(value).lower() == "true":
        return True
    elif str(value).lower() == "false":
        return False
    else:
        raise ValueError(f"cannot cast {value} as boolean")


def value_is_nan(value):
    """Check if value is nan
    nan is a special dataframe value to indicate missing data
    """
    # We coerces nan to string for group annotations (SCP-2545)
    # string value 'nan' should also be considered value_is_nan=True
    if str(value).lower() == "n/a":
        return True
    try:
        number = float(value)
        return math.isnan(number)
    except TypeError:
        return False
    except ValueError:
        return False


def cast_integer_type(value):
    """Cast metadata value as integer"""
    if value_is_nan(value):
        # nan indicates missing data, has no valid integer value for casting
        return value
    else:
        return int(value)


def cast_float_type(value):
    """Cast metadata value as float"""
    return float(value)


def cast_string_type(value):
    """Cast string type per convention where Pandas autodetected a number"""
    if value_is_nan(value):
        # nan indicates missing data; by type, nan is a numpy float
        # so a separate type check is needed for proper handling
        return value
    elif isinstance(value, numbers.Number):
        return str(value)
    else:
        return value


def regularize_ontology_id(value):
    """Regularize ontology_ids for storage with underscore format"""
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
        "number": cast_float_type,
        "boolean": cast_boolean_type,
        "integer": cast_integer_type,
        "string": cast_string_type,
    }
    if is_array_metadata(convention, metadatum):
        cast_values = []
        # check
        metadata.update_numeric_array_columns(metadatum)
        try:
            if "|" not in value:
                msg = (
                    f"{metadatum}: accepts an array of values.\n"
                    "Lines detected with single values instead of arrays.\n"
                    "If multiple values are expected, use a pipe ('|') to separate values."
                )
                metadata.store_validation_issue(
                    "warn",
                    msg,
                    "content:array-no-pipes",
                    associated_info=[id_for_error_detail],
                )

            # splitting on pipe character for array data, valid for Sarah's
            # programmatically generated SCP TSV metadata files. When ingesting
            # files that support array-based metadata navtively (eg. loom,
            # anndata etc) splitting on pipe may become problematic
            for element in value.split("|"):
                try:
                    if "ontology" in convention["properties"][metadatum]:
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
            msg = (
                f"{metadatum}: '{element}' in '{value}' does not match "
                f"expected '{lookup_metadata_type(convention, metadatum)}' type."
            )
            metadata.store_validation_issue(
                "error",
                msg,
                "content:type:value-type-mismatch",
                associated_info=[id_for_error_detail],
            )
        # This exception should only trigger if a single-value array
        # metadata is being cast - the value needs to be passed as an array,
        except (AttributeError, TypeError):
            try:
                if "ontology" in convention["properties"][metadatum]:
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
            if "ontology" in convention["properties"][metadatum]:
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
            msg = f'{metadatum}: "{value}" does not match expected type.'
            metadata.store_validation_issue(
                "error",
                msg,
                "content:type:value-type-mismatch",
                associated_info=[id_for_error_detail],
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
    metadata_names[0] = "CellID"
    row_info = dict(itertools.zip_longest(metadata_names, line))
    processed_row = {}
    for k, v in row_info.items():
        # for metadata not in convention, check numeric metadata for valid values
        if k not in convention["properties"].keys():
            type_index = metadata.headers.index(k)
            k_type = metadata.annot_types[type_index]
            if k_type == "numeric":
                k_numeric = type(v) in [int, float]
                if not k_numeric:
                    # pandas coercion will stringify numbers
                    # error messages including v are more confusing than helpful
                    # e.g. percent_mt: supplied value 0.0825991189427313 is not numeric
                    # 0.0825991189427313 was stringified due to non-numeric in the data column
                    # but the value 0.0825991189427313 was validly supplied as a numeric
                    try:
                        float(v)
                        msg = f"{k}: one or more values in data column are not numeric."
                    # only true non-numerics should be reported in detail
                    except ValueError:
                        msg = f"{k}: supplied value {v} is not numeric."
                    metadata.store_validation_issue(
                        "error",
                        msg,
                        "content:type:not-numeric",
                        associated_info=row_info["CellID"],
                    )
                    dev_logger.error(msg)
            continue
        # for optional metadata, do not pass empty cells (nan)
        if k not in convention["required"]:
            if value_is_nan(v) or is_empty_string(v):
                continue
        else:
            if is_array_metadata(convention, k):
                for element in v.split("|"):
                    if value_is_nan(element) or is_empty_string(element):
                        msg = f"{k}: NA, NaN, and None are not accepted values for arrays or required fields."
                        metadata.store_validation_issue(
                            "error",
                            msg,
                            "content:missing-required-values",
                            associated_info=row_info["CellID"],
                        )
                        dev_logger.error(msg)

        processed_row.update(
            cast_metadata_type(k, v, row_info["CellID"], convention, metadata)
        )
    return processed_row


def is_empty_string(value):
    if isinstance(value, str):
        return value == "" or value.isspace()


def collect_jsonschema_errors(metadata, convention, bq_json=None):
    """Evaluate metadata input against metadata convention using JSON schema
    returns False if input convention is invalid JSON schema
    """
    # this function seems overloaded with its three tasks
    # schema validation, non-ontology errors, ontology info collection
    # the latter two should be done together in the same pass through the file
    js_errors = defaultdict(list)
    schema = validate_schema(convention, metadata)

    if bq_json:
        bq_filename = str(metadata.study_file_id) + ".json"
        # truncate JSON file so data from serialize_bq starts with an empty file
        fh = open(bq_filename, "w")
        fh.close()
    if schema:
        compare_type_annots_to_convention(metadata, convention)
        rows = metadata.yield_by_row()
        line = next(rows)
        while line:
            row = process_metadata_row(metadata, convention, line)
            metadata.cells.append(row["CellID"])
            row = collect_ontology_data(row, metadata, convention)
            for error in schema.iter_errors(row):
                try:
                    error.message = error.path[0] + ": " + error.message
                except IndexError:
                    pass
                js_errors[error.message].append(row["CellID"])
            if bq_json:
                # add non-convention, SCP-required, metadata for BigQuery
                row["study_accession"] = metadata.study_accession
                row["file_id"] = str(metadata.study_file_id)
                # extract convention version from URI in convention JSON file
                row["metadata_convention_version"] = convention["$id"].split("/")[-2]
                serialize_bq(row, bq_filename)
            try:
                line = next(rows)
            except StopIteration:
                break
        if js_errors:
            metadata.store_validation_issue(
                "error",
                "One or more errors detected validating metadata content against metadata convention",
                "convention:jsonschema-error",
            )
            metadata.issues["error"]["convention"] = js_errors
        validate_cells_unique(metadata)
    else:
        return False


def record_issue(issue_type, msg):
    """Print issue to console with coloring and
    send errors to logger
    """

    if issue_type == "error":
        user_logger.error(msg)
        color = Fore.RED
    elif issue_type == "warn":
        color = Fore.YELLOW
    else:
        color = ""
    console_msg = color + msg
    print(console_msg)


def report_issues(metadata):
    """Report issues in CellMetadata.issues dictionary
    returns True if errors are reported, False if no errors to report
    """

    dev_logger.info(f"Checking for validation issues")

    has_errors = False
    has_warnings = False
    for issue_type in sorted(metadata.issues.keys()):
        for issue_category, category_dict in metadata.issues[issue_type].items():
            if category_dict:
                category_header = f"\n*** {issue_category} {issue_type} list:"
                record_issue(issue_type, category_header)
                if issue_type == "error":
                    has_errors = True
                elif issue_type == "warn":
                    has_warnings = True
                for issue_text, cells in category_dict.items():
                    if cells:
                        try:
                            re.compile(issue_text)
                            is_valid = True
                        except re.error:
                            pass
                        if issue_type == "error":
                            issue_msg = f"{issue_text} [ Error count: {len(cells)} ]"
                        elif issue_type == "warn":
                            issue_msg = f"{issue_text} [ Warn count: {len(cells)} ]"
                        record_issue(issue_type, issue_msg)
                    else:
                        record_issue(issue_type, issue_text)
    if not has_errors and not has_warnings:
        no_issues = "No errors or warnings detected for input metadata file"
        record_issue(None, no_issues)
    return has_errors


def exit_if_errors(metadata):
    """Determine if CellMetadata.issues has errors
    Exit with error code 1 if errors are reported, return False if no errors
    """
    errors = False
    for error_type in metadata.issues.keys():
        for error_category, category_dict in metadata.issues[error_type].items():
            if category_dict:
                if error_type == "error":
                    errors = True
    if errors:
        exit(1)
    return errors


def is_label_or_synonym(labels, provided_label):
    """Determine if a user-provided ontology label is a valid label or synonymn
    :param labels: cached ontology label/synonyms from retriever.retrieve_ontology_term_label_and_synonyms
    :param provided_label: user-provided label from metadata file
    :return: True/False on match for label or synonym
    """
    label = labels.get('label')
    synonyms = labels.get('synonyms')
    if label.casefold() == provided_label.casefold():
        return True
    elif synonyms:
        if next((t for t in synonyms if t.casefold() == provided_label.casefold()), ''):
            return True
        else:
            return False
    else:
        return False

def ontology_is_local(ontology_name):
    """Check if it is possible to use local ontology validation instead of OLS
    :param ontology_name: name of ontology
    :return: Boolean
    """
    return ontology_name is not None and ontology_name in minified_reader.ontology_names()

def validate_collected_ontology_data(metadata, convention):
    """Evaluate collected ontology_id, ontology_label info in
    CellMetadata.ontology dictionary by querying ontology source
    validity of ontology_id and cross-check that input ontology_label
    matches ontology_label from ontology lookup
    """
    for property_name in metadata.ontology.keys():
        # pull in file for validation of non-EBI OLS metadata for organ_region
        if property_name == "organ_region":
            try:
                region_ontology = (  # noqa: F841
                    retriever.fetch_allen_mouse_brain_atlas()
                )
            except BaseException as e:
                print(e)
                # immediately return as validation should not continue if ontology file unavailable
                return
        # split on comma in case this property from the convention supports multiple ontologies

        ontology_urls = convention["properties"][property_name]["ontology"].split(",")

        for ontology_info in metadata.ontology[property_name].keys():
            ontology_id, ontology_label = ontology_info
            try:
                attribute_type = convention["properties"][property_name]["type"]
                # get actual label along with synonyms for more robust matching
                label_and_synonyms = retrieve_label_and_synonyms(ontology_id, property_name, convention, attribute_type)

                if not is_label_or_synonym(label_and_synonyms, ontology_label):
                    matched_label_for_id = label_and_synonyms.get("label")
                    ontology_source_name = "EBI OLS"
                    if property_name == "organ_region":
                        ontology_source_name = "Allen Mouse Brain Atlas"
                    if is_required_metadata(convention, property_name) and (
                        not ontology_label or value_is_nan(ontology_label)
                    ):
                        msg = (
                            f'{property_name}: ontology label required for "{ontology_id}" '
                            f'to cross-check for data entry error.'
                        )
                        metadata.store_validation_issue(
                            "error",
                            msg,
                            "content:missing-required-values",
                            associated_info=metadata.ontology[property_name][
                                (ontology_id, ontology_label)
                            ],
                        )
                    else:
                        msg = (
                            f'{property_name}: input ontology_label "{ontology_label}" '
                            f'does not match {ontology_source_name} lookup "{matched_label_for_id}" for ontology id "{ontology_id}".'
                        )
                        metadata.store_validation_issue(
                            "error",
                            msg,
                            "ontology:label-not-match-id",
                            associated_info=metadata.ontology[property_name][
                                (ontology_id, ontology_label)
                            ],
                        )
                else:
                    property_header = property_name + "__ontology_label"
                    matched_label_for_id = label_and_synonyms.get("label")
                    if ontology_label != matched_label_for_id:
                        metadata.synonym_updates[property_header][
                            ontology_label
                        ] = matched_label_for_id
            except ValueError as value_error:
                metadata.store_validation_issue(
                    "error",
                    value_error.args[0],
                    "ontology:label-lookup-error",
                    associated_info=metadata.ontology[property_name][
                        (ontology_id, ontology_label)
                    ],
                )
            except requests.exceptions.RequestException as err:
                msg = f"External service outage connecting to {ontology_urls} when querying {ontology_id}:{ontology_label}: {err}"
                dev_logger.exception(msg)
                metadata.store_validation_issue(
                    "error", msg, "runtime:request-exception", issue_type="runtime"
                ),
                # immediately return as validation cannot continue
                return None

    return


def confirm_uniform_units(metadata, convention):
    """Check that any unit metadata are uniform within study
    Note: refactoring may be needed if metadata files are chunked
    """
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    for name in metadata_names:
        if name.endswith("__unit"):
            # existence of unit metadata is enforced by jsonschema
            # any empty cells in a unit column are okay to be empty
            if metadata.file[name].nunique(dropna=True).values[0] != 1:
                msg = f"{name}: values for each unit metadata required to be uniform."
                metadata.store_validation_issue("error", msg, "content:uniform-units")


def serialize_bq(bq_dict, filename="bq.json"):
    """Write metadata collected for validation to json file
    BigQuery requires newline delimited json objects
    """
    bq_dict_copy = bq_dict.copy()
    if "organism_age" in bq_dict and "organism_age__unit_label" in bq_dict:
        bq_dict_copy["organism_age__seconds"] = calculate_organism_age_in_seconds(
            bq_dict["organism_age"], bq_dict["organism_age__unit_label"]
        )

    data = json.dumps(bq_dict_copy)
    with open(filename, "a") as jsonfile:
        jsonfile.write(data + "\n")


def calculate_organism_age_in_seconds(organism_age, organism_age_unit_label):
    multipliers = {
        "microsecond": 0.000001,
        "millisecond": 0.001,
        "second": 1,
        "minute": 60,
        "hour": 3600,
        "day": 86400,
        "week": 604800,
        "month": 2626560,  # (day * 30.4 to fuzzy-account for different months)
        "year": 31557600,  # (day * 365.25 to fuzzy-account for leap-years)
    }
    multiplier = multipliers[organism_age_unit_label]
    return organism_age * multiplier


def serialize_issues(metadata):
    """Write collected issues to json file"""
    with open("issues.json", "w") as jsonfile:
        json.dump(metadata.issues, jsonfile, indent=2)
        jsonfile.write("\n")  # Add newline cause Py JSON does not


def review_metadata_names(metadata):
    """Check metadata names for disallowed characters"""
    metadata_names = metadata.file.columns.get_level_values(0).tolist()
    for name in metadata_names:
        allowed_char = re.compile("^[A-Za-z0-9_]+$")
        if not allowed_char.match(name):
            msg = (
                f"{name}: only alphanumeric characters and underscore "
                f"allowed in metadata name."
            )
            metadata.store_validation_issue(
                "error", msg, "format:cap:only-alphanumeric-underscore"
            )


def identify_multiply_assigned(list):
    """Given a list of ontology IDs and their purported labels,
    return list of unique multiply-assigned labels
    """
    ontology_tracker = defaultdict(lambda: defaultdict(int))
    multiply_assigned = []
    for element in list:
        id, label = element
        ontology_tracker[label][id] += 1
    for label in ontology_tracker:
        if len(ontology_tracker[label].keys()) > 1:
            multiply_assigned.append(label)
    multiply_assigned.sort()
    return multiply_assigned


def assess_ontology_ids(ids, property_name, metadata):
    """
    Check ordered collection of ontology IDs for increasing numeric values
    """
    evidence_of_excel_drag = False
    evidence_of_excel_drag_threshold = 25
    binned_ids = defaultdict(list)
    for id in ids:
        # The binning avoids any spurrious numeric contiguity
        # between IDs that actually have different shortnames
        # because the detection threshold is fairly generous,
        # the binning is probably unneeded
        try:
            ontology_shortname, term_id = re.split("[_:]", id)
            binned_ids[ontology_shortname].append(term_id)
        except (ValueError, TypeError):
            # invalid ontology id will already be flagged as convention error
            # no need to also mark as ontology error as part of Excel drag detection
            pass
    for ontology in binned_ids.keys():
        id_numerics = []
        incrementation_count = 0
        for term in binned_ids[ontology]:
            # Regex extracts away the numeric part of the term
            # from the constant, text portion of the ontology ID
            # some term ids have a text component on the term side
            # so a regex is needed instead of a simple split
            term_numeric = re.search('(\d)*$', term)
            id_numerics.append(int(term_numeric.group()))
        for x, y in zip(id_numerics, id_numerics[1:]):
            if y - x == 1:
                incrementation_count += 1
            elif y - x != 1:
                incrementation_count = 0
            if incrementation_count >= evidence_of_excel_drag_threshold:
                evidence_of_excel_drag = True
                break
    return evidence_of_excel_drag


def detect_excel_drag(metadata, convention):
    """Check if ontology IDs submitted have characteristic Excel drag properties
    Todo1: "Excel drag" detection of array-based ontologyID data
    is lacking (need to track pipe-delimited string)
    Todo2: need to bypass EBI OLS queries to "fill in"
    missing ontology labels for optional metadata)
    Hint: try working with raw metadata.file data, would
    allow this check to be moved before collect_jsonschema_errors
    """
    excel_drag = False
    for property_name in metadata.ontology.keys():
        if len(set(metadata.ordered_ontology[property_name])) == 1:
            continue
        else:
            property_ids = metadata.ordered_ontology[property_name]
            unique_ids = set(property_ids)
            property_labels = metadata.ordered_labels[property_name]
            property_labels_blanks_removed = [i for i in property_labels if i]
            unique_labels = set(property_labels_blanks_removed)
            # likely ontology label mis-assignment if multiple ontology IDs ascribed to same ontology label
            label_multiply_assigned = len(unique_ids) > len(unique_labels)
            try:
                if assess_ontology_ids(property_ids, property_name, metadata):
                    msg = (
                        f"{property_name}: Long stretch of contiguously incrementing "
                        + "ontology ID values suggest cut and paste issue - exiting validation, "
                        + "ontology content not validated against ontology server.\n"
                        "Please confirm ontology IDs are correct and resubmit.\n"
                    )
                    excel_drag = True
                    if label_multiply_assigned:
                        multiply_assigned = identify_multiply_assigned(
                            set(zip(property_ids, property_labels))
                        )
                        if multiply_assigned:
                            msg += f"Check for mismatches between ontology ID and provided ontology label(s) {multiply_assigned}.\n"
                    metadata.store_validation_issue(
                        "error", msg, "ontology:multiply-assigned-label"
                    )
                    dev_logger.exception(msg)
            except ValueError as value_error:
                metadata.store_validation_issue(
                    "error", value_error.args[0], "ontology:label-lookup-error"
                )

    return excel_drag


def to_1D(series):
    """Pandas values which are list need unnesting"""
    return [x for _list in series for x in _list]


def replace_single_value_array(df, metadata_name, synonym, label):
    """Synonym replacement (in-place) for single-value array metadata
    Pandas doesn't operate well on lists which are potentially non-homogenous
    # https://stackoverflow.com/questions/53116286/how-to-assign-an-entire-list-to-each-row-of-a-pandas-dataframe
    """
    match = [v == [synonym] for v in df[metadata_name]]
    value = [label]
    df.loc[match, metadata_name] = df.apply(lambda x: value, axis=1)


def replace_synonym_in_multivalue_array(df, metadata_name, substitutions):
    """Synonym replacement (in-place) for multi-value array ontology labels
    must identify all affected arrays of labels, construct replacement arrays
    then replace old synonym-containing array with an updated array
    """
    orig_values = list(df[metadata_name].transform(tuple).unique())
    matching_synonyms = {}
    for o in orig_values:
        # if a synonym (s) is an element of the multivalue array (o), track it
        matching_synonyms[o] = [s for s in substitutions.keys() if s in o]
    for o in matching_synonyms.keys():
        # if a multivalue array contains synonyms
        if matching_synonyms[o]:
            # make a copy of the original multivalue array that will take on all substitutions
            replacement_value = list(o)
            # make all synonym substitutions into the multivalue array of ontology labels
            for s in matching_synonyms[o]:
                replacement_value = [
                    substitutions[s] if term == s else term
                    for term in replacement_value
                ]
            # select the rows in the dataframe with entries for multivalue array (o)
            match = [v == list(o) for v in df[metadata_name]]
            df.loc[match, metadata_name] = df.apply(lambda x: replacement_value, axis=1)


def replace_synonyms(metadata):
    """
    Update BigQuery data to store ontology labels and not synonyms
    """
    bq_filename = str(metadata.study_file_id) + ".json"
    df = pd.read_json(bq_filename, lines=True)
    for metadata_name in metadata.synonym_updates.keys():
        # non-array metadata values are strings
        if isinstance(df[metadata_name][0], str):
            for synonym in metadata.synonym_updates[metadata_name].keys():
                df[metadata_name].replace(
                    synonym,
                    metadata.synonym_updates[metadata_name][synonym],
                    inplace=True,
                )
        # Pandas can't hash mutable complex objects like lists
        # need to find and replace by location (iloc)
        elif len(df[metadata_name]) == len(to_1D(df[metadata_name])):
            for synonym in metadata.synonym_updates[metadata_name].keys():
                replace_single_value_array(
                    df,
                    metadata_name,
                    synonym,
                    metadata.synonym_updates[metadata_name][synonym],
                )
        # at least one non-single array-type metadata
        else:
            replace_synonym_in_multivalue_array(
                df, metadata_name, metadata.synonym_updates[metadata_name]
            )
    df.to_json(bq_filename, orient="records", lines=True)


def validate_input_metadata(metadata, convention, bq_json=None):
    """Wrapper function to run validation functions"""
    dev_logger.info("Checking metadata content against convention rules")
    collect_jsonschema_errors(metadata, convention, bq_json)
    review_metadata_names(metadata)
    dev_logger.info('Checking for "Excel drag" events')
    if not detect_excel_drag(metadata, convention):
        # "short-circuit" ontology validation if "Excel drag" detected
        # avoids a bloat of calls to EBI OLS, return error faster and avoid
        # long-compute-time issue (if false positives are possible, bypass will be needed)
        dev_logger.info('Validating ontology content against EBI OLS')
        validate_collected_ontology_data(metadata, convention)
        if metadata.synonym_updates:
            replace_synonyms(metadata)
        confirm_uniform_units(metadata, convention)


def push_metadata_to_bq(metadata, ndjson, dataset, table):
    """Upload local NDJSON to BigQuery"""
    client = bigquery.Client()
    dataset_ref = client.dataset(dataset)
    table_ref = dataset_ref.table(table)
    job_config = bigquery.LoadJobConfig()
    job_config.write_disposition = bigquery.WriteDisposition.WRITE_APPEND
    job_config.source_format = bigquery.SourceFormat.NEWLINE_DELIMITED_JSON
    try:
        with open(ndjson, "rb") as source_file:
            job = client.load_table_from_file(
                source_file, table_ref, job_config=job_config
            )
        job.result()  # Waits for table load to complete.
        info_msg = f"Metadata uploaded to BigQuery. ({job.output_rows} rows)"
        print(info_msg)
        dev_logger.info(info_msg)
    # Unable to intentionally trigger a failed BigQuery upload
    # please add print statement below to error logging so
    # error handling can be updated when better understood
    except Exception as e:
        print(e)
        dev_logger.exception(e)
        return 1
    if job.output_rows != len(metadata.cells):
        msg = f"BigQuery upload error: upload ({job.output_rows} rows) does not match number of cells in file, {len(metadata.cells)} cells."
        print(msg)
        dev_logger.error(msg)
        return 1
    os.remove(ndjson)
    return 0


def write_metadata_to_bq(metadata, bq_dataset, bq_table):
    """Wrapper function to gather metadata and write to BigQuery"""
    bq_filename = str(metadata.study_file_id) + ".json"
    push_status = push_metadata_to_bq(metadata, bq_filename, bq_dataset, bq_table)
    return push_status


def check_if_old_output():
    """Exit if old output files found"""
    output_files = ["bq.json"]

    old_output = False
    for file in output_files:
        if os.path.exists(file):
            print(f"{file} already exists, please delete file and try again")
            old_output = True
    if old_output:
        exit(1)
