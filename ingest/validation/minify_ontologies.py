"""Minifies ontologies used in EBI OLS, to enable instant ontology validation

This converts ~224 MB in ontology JSON files into 2 MB TSV.GZs at build-time.
The 2 MB compressed ontologies can then be retrieved at runtime.
Only IDs, labels, and synonyms are retained from the original ontologies.

Example:
cd ingest/validation
python minify_ontologies.py
"""

import argparse
import json
import urllib.request
from pathlib import Path
import gzip

MONDO_URL = 'https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.json'
PATO_URL = 'https://github.com/pato-ontology/pato/raw/master/pato.json'
NCBITAXON_URL = 'https://github.com/obophenotype/ncbitaxon/releases/latest/download/taxslim.json'
EFO_URL = 'https://github.com/EBISPOT/efo/releases/latest/download/efo.json'
UBERON_URL = 'https://github.com/obophenotype/uberon/releases/latest/download/uberon.json'

ONTOLOGY_JSON_URLS = {
    'disease': [MONDO_URL, PATO_URL],
    'species': [NCBITAXON_URL],
    'library_preparation_protocol': [EFO_URL],
    'organ': [UBERON_URL]
}

def fetch(url, use_cache=True):
    """Request remote resource, read local cache if availalble
    """
    filename = url.split('/')[-1]
    if use_cache == False or (use_cache and not Path(filename).is_file()):
        with urllib.request.urlopen(url) as f:
            content = f.read()
        if use_cache:
            with open(filename, 'wb') as f:
                f.write(content)
    else:
        with open(filename) as f:
            content = f.read()
    return [content, filename]

def fetch_ontologies(ontology_json_urls, use_cache=True):
    """Retrieve ontology JSON and JSON filename for required ontology
    """
    ontologies = {}
    for annotation in ontology_json_urls:
        ontology_urls = ontology_json_urls[annotation]
        ontologies[annotation] = []
        for ontology_url in ontology_urls:
            print(f'Fetch ontology: {ontology_url}')
            raw_ontology, filename = fetch(ontology_url, use_cache)
            ontology_json = json.loads(raw_ontology)
            ontologies[annotation].append([ontology_json, filename])
    return ontologies

def get_synonyms(node, label):
    """Get related and exact synonyms for an ontology node
    """
    if 'meta' not in node or 'synonyms' not in node['meta']:
        return ''

    raw_synonyms = []
    synonym_nodes = node['meta']['synonyms']
    for synonym_node in synonym_nodes:
        if 'val' not in synonym_node:
            # Handles e.g. incomplete EFO synonym nodes
            continue
        raw_synonym = synonym_node['val']
        if (
            not raw_synonym.startswith('obsolete ') and # Omit obsolete synonyms
            raw_synonym != label # Omit synonyms that are redundant with label
        ):
            raw_synonyms.append(raw_synonym)
    synonyms = '||'.join(raw_synonyms) # Unambiguously delimit synonyms
    return synonyms

def minify(ontology_json, filename):
    """Convert full ontology JSON into a minimal gzipped TSV, write to disk
    """
    ontology_shortname = filename.split('.json')[0]
    if ontology_shortname == 'taxslim':
        ontology_shortname = 'ncbitaxon'
    ontology_shortname_uc = ontology_shortname.upper()
    graph_nodes = ontology_json['graphs'][0]['nodes']

    raw_nodes = list(filter(
        lambda n: f'/{ontology_shortname_uc}_' in n['id'].upper() and 'lbl' in n,
        graph_nodes
    ))

    all_nodes = list(map(
        lambda n: (
            [n['id'].split('/')[-1], n['lbl'], get_synonyms(n, n['lbl'])]
        ), raw_nodes
    ))

    # Remove obsolete labels
    nodes = list(filter(
        lambda n: not n[1].startswith('obsolete '),
        all_nodes
    ))

    tsv_content = '\n'.join(
        map(lambda n: '\t'.join(n), nodes)
    )
    compressed_tsv_content = gzip.compress(tsv_content.encode())

    output_filename = f'ontologies/{ontology_shortname}.min.tsv.gz'
    with open(output_filename, 'wb') as f:
        f.write(compressed_tsv_content)
    print(f'Wrote {output_filename}')


class OntologyMinifier:

    def __init__(self, annotations=None, use_cache=True):
        # Enable minifying incomplete set of ontologies, e.g. for testing
        if annotations:
            ontology_json_urls = {}
            for annotation in annotations:
                ontology_json_urls[annotation] = ONTOLOGY_JSON_URLS[annotation]
        else:
            ontology_json_urls = ONTOLOGY_JSON_URLS

        ontologies = fetch_ontologies(ontology_json_urls, use_cache)
        for annotation in ontologies:
            for conf in ontologies[annotation]:
                ontology_json, filename = conf
                minify(ontology_json, filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--use-cache",
        help=(
            "Whether to use previously-downloaded raw ontologies"
        ),
        action="store_true"
    )
    args = parser.parse_args()
    use_cache = args.use_cache
    OntologyMinifier(None, use_cache)
