"""Command-line interface for validating ontology-based terms/IDs against EBI Ontology Lookup Service

DESCRIPTION
This CLI will read a tab- or comma-delimited cell metadata file and validate any columns that are
designated as ontologically based against the EBI OLS.  Can also retrieve all ancestors for a given term.

PREREQUISITES
Requests 2.22.0 (see requirements.txt)

EXAMPLES
Validate columns called 'Cell Type' and 'Disease' from a cell metadata file against EBI OLS
python3 validate_ontology_terms.py --metadata-path /path/to/metadata.txt --columns "Cell Type, Disease"
"""

import argparse
import time
import os
import requests
import urllib.parse as encoder

start_time = time.time()

OLS_BASE_URL = "https://www.ebi.ac.uk/ols/api/ontologies/"
terms_not_found = {}
error_filename = "metadata_errors.txt"

parser = argparse.ArgumentParser(
    prog='validate_ontology_terms.py'
)

parser.add_argument(
    '-m', '--metadata-path', default=None,
    help='Path to tsv/csv cell metadata file'
)

parser.add_argument(
    '-c', '--columns', default=None,
    help='Comma-delimited list of specific columns to validate from cell metadata file'
)

parser.add_argument(
    '-i', '--ignore-columns', default="Cell ID",
    help='Comma-delimited list of specific columns to ignore from cell metadata file'
)

parser.add_argument(
    '-d', '--delimiter', default="\t",
    help='Delimiter for cell metadata file ("\\t" as default)'
)

parser.add_argument(
    '-a', '--ancestors', default=False, action='store_true',
    help='Retrieve ancestors for ontology terms from defining ontology (False as default)'
)

args = parser.parse_args()

# ensure cell metadata file exists
if args.metadata_path:
    metadata_exists = os.path.isfile(args.metadata_path)
    if not metadata_exists:
        print("Metadata file not found at: " + args.metadata_path)
        exit(1)
else:
    print("No metadata file provided; Exiting")
    exit(1)


# set control lists of columns to validate/ignore
validation_columns = []
excluded_columns = []
if args.columns:
    controls = args.columns.split(',')
    validation_columns = [column.strip() for column in controls]
    print("Validating declared columns: %s." % validation_columns)

if args.ignore_columns:
    ignores = args.ignore_columns.split(',')
    excluded_columns = [column.strip() for column in ignores]
    print("Excluding columns: %s." % excluded_columns)


def load_cell_metadata(filepath, control_list, exclusion_list):
    """Read cell metadata from a plain text file and return dict of unique values by column.
    :param filepath: filepath to input cell metadata file
    :param control_list: list of columns to validate (exclusive - will ignore all other columns)
    :param exclusion_list: list of columns to exclude
    :return: dict of column names and unique terms
    """
    metadata = {}
    with open(args.metadata_path) as file:
        headers = file.readline().rstrip("\n").split(args.delimiter)
        for header in headers:
            metadata[header] = []
        # remove ignored columns
        for ignore in excluded_columns:
            if ignore in metadata:
                del metadata[ignore]
        # only validate declared columns
        if validation_columns:
            cols_to_remove = list(set(metadata.keys()) - set(validation_columns))
            for col in cols_to_remove:
                del metadata[col]
        for line in file:
            values = line.rstrip("\n").split(args.delimiter)
            for idx, header in enumerate(headers):
                value = values[idx]
                if header in metadata and value not in metadata[header]:
                    metadata[header].append(value)
    return metadata


def extract_identifiers_from_term(ontology_term):
    """ Extract the ontology identifier and term ID from a term (e.g. CL_0002419 returns "CL" and "0002419")
    :param ontology_term: ontology term id containing both an ontology ID and term ID delimited by _
    :return: list of ontology ID and term ID
    """
    return ontology_term.split('_')


def retrieve_ontology(ontology_term):
    """Retrieve an ontology listing from EBI OLS
    :param ontology_term: identifier of a term in an ontology in OLS (e.g. CL_0002419)
    :return: JSON payload of ontology, or None
    """
    ontology_id, term_id = extract_identifiers_from_term(ontology_term)
    response = requests.get(OLS_BASE_URL + ontology_id)
    if response.status_code == 200:
        return response.json()
    else:
        return None


def retrieve_ontology_term(ontology_term):
    """Retrieve an individual term from an ontology
    :param ontology_term: term to query for in matching ontology
    :return: JSON payload of ontology of ontology term, or None
    """
    ontology = retrieve_ontology(ontology_term)
    if ontology:
        base_term_uri = ontology['config']['baseUris'][0]
        query_iri = encode_term_iri(ontology_term, base_term_uri)
        term_url = ontology['_links']['terms']['href'] + '/' + query_iri
        response = requests.get(term_url)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    else:
        return None


def retrieve_ontology_term_ancestors(ontology_term):
    """Retrieve all ancestors for a given term, scoped to definiting ontology

    :param ontology_term: identifier of a term in an ontology in OLS (e.g. CL_0002419)
    :return: list of ontology IDs for all direct ancestors
    """
    term_data = retrieve_ontology_term(ontology_term)
    if term_data:
        ancestors_url = term_data['_links']['ancestors']['href']
        response = requests.get(ancestors_url)
        if response.status_code == 200:
            ancestors = []
            all_ancestors = response.json()
            for ancestor in all_ancestors['_embedded']['terms']:
                if ancestor['is_defining_ontology']:
                    ancestors.append(ancestor['iri'].split('/')[-1])
            return ancestors
        else:
            return None


def encode_term_iri(term, base_uri):
    """Double url-encode a term Internationalized Resource Identifier (IRI) for querying OLS ontologies

    :param term: ontology term
    :param base_uri: base term URI for corresponding ontology
    :return: double url-encoded ontology term IRI
    """
    ontology_id, term_id = extract_identifiers_from_term(term)
    ref_ontology_id = base_uri.split('/')[-1]
    if ontology_id == ref_ontology_id:
        query_uri = base_uri + term_id
    else:
        query_uri = base_uri.replace("/" + ref_ontology_id, "/") + term
    encoded_iri = encoder.quote_plus(encoder.quote_plus(query_uri))
    return encoded_iri


def output_error_report(invalid_terms):
    """Output a list of missing terms, organized by column

    :param invalid_terms: list of terms not found in OLS
    """
    with open(error_filename, 'w') as error_file:
        error_file.write("Column Name%sInvalid Terms\n" % args.delimiter)
        for row_name, values in invalid_terms.items():
            row = row_name + args.delimiter
            row += args.delimiter.join(values)
            error_file.write(row + "\n")


print("Ingesting cell metadata...")
cell_metadata = load_cell_metadata(args.metadata_path, validation_columns, excluded_columns)
print("Ingest complete! Detected the following unique values")
for column in cell_metadata.keys():
    print("Column: " + column)
    print(cell_metadata[column])

total_matches = 0
print("Validating terms against OLS...")
for column, term_list in cell_metadata.items():
    for term in term_list:
        print("####")
        print("")
        print("Checking for existence of " + term)
        matching_term = retrieve_ontology_term(term)
        if matching_term:
            total_matches += 1
            print("Found matching for " + term)
            print("")
            print("### DATA FOR " + term + " ###")
            print("")
            print("Label: " + matching_term['label'])
            if matching_term['description']:
                print("Description: " + matching_term['description'][0])
            if args.ancestors:
                print('Retrieving ancestors...')
                term_ancestors = retrieve_ontology_term_ancestors(term)
                if term_ancestors:
                    print("Ancestors: %s" % ', '.join(term_ancestors))
        else:
            print("No match found for " + term)
            if column not in terms_not_found.keys():
                terms_not_found[column] = []
            if term not in terms_not_found[column]:
                terms_not_found[column].append(term)
        print("")
        print("####")

total_terms = sum(len(vals) for vals in cell_metadata.values())
total_missing = sum(len(vals) for vals in terms_not_found.values())
print("Validation complete! Matched %s of %s terms" % (total_matches, total_terms))
total_time = str(round(time.time() - start_time))
print('Total time: ' + total_time + 's')
if terms_not_found:
    print("%s terms not found, writing errors to %s" % (total_missing, error_filename))
    output_error_report(terms_not_found)
    exit(1)
else:
    exit(0)


