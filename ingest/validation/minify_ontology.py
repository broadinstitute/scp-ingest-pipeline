import json
import urllib
import urllib.request
from pathlib import Path

mondo_url = 'https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.json'
pato_url = 'https://github.com/pato-ontology/pato/raw/master/pato.json'

# TODO: Determine how to query for current JSON
ncbitaxon_url = 'https://github.com/obophenotype/ncbitaxon/releases/latest/download/taxslim.json'
efo_url = 'https://github.com/EBISPOT/efo/releases/latest/download/efo.json'

ONTOLOGY_JSON_URLS = {
    'disease': [mondo_url, pato_url],
    'species': [ncbitaxon_url],
    'library_preparation_protocol': [efo_url]
}

def fetch(url, use_cache=True):
    """Request remote resource, read local cache if availalble
    """
    filename = url.split('/')[-1]
    if use_cache and not Path(filename).is_file():
        with urllib.request.urlopen(url) as f:
            content = f.read()
        with open(filename, 'wb') as f:
            f.write(content)
    else:
        with open(filename) as f:
            content = f.read()
    return [content, filename]


def fetch_ontologies(use_cache=True):
    """Retrieve ontology JSON and JSON filename for required ontology
    """
    ontologies = {}
    for annotation in ONTOLOGY_JSON_URLS:
        ontology_urls = ONTOLOGY_JSON_URLS[annotation]
        ontologies[annotation] = []
        for ontology_url in ontology_urls:
            print(f'Fetch ontology: {ontology_url}')
            raw_ontology, filename = fetch(ontology_url, use_cache)
            ontology_json = json.loads(raw_ontology)
            ontologies[annotation].append([ontology_json, filename])
    return ontologies


def get_synonyms(node):
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
        raw_synonyms.append(synonym_node['val'])
    # print('raw_synonyms', raw_synonyms)
    synonyms = '||'.join(raw_synonyms) # Unambiguously delimit synonyms
    return synonyms

def minify(ontology_json, filename):
    print(f'Minify {filename}')
    ontology_shortname = filename.split('.json')[0]
    if ontology_shortname == 'taxslim':
        ontology_shortname = 'ncbitaxon'
    ontology_shortname_uc = ontology_shortname.upper()
    graph_nodes = ontology_json['graphs'][0]['nodes']


    raw_nodes = list(filter(
        lambda n: f'/{ontology_shortname_uc}_' in n['id'].upper() and 'lbl' in n,
        graph_nodes
    ))

    nodes = list(map(
        lambda n: '\t'.join(
            [n['id'].split('/')[-1], n['lbl'], get_synonyms(n)]
        ), raw_nodes
    ))

    with open(f'{ontology_shortname}.min.tsv', 'w') as f:
        f.write('\n'.join(nodes))

def run(use_cache=True):
    print('Run')
    ontologies = fetch_ontologies(use_cache)
    for annotation in ontologies:
        for conf in ontologies[annotation]:
            # print('conf', conf)
            ontology_json, filename = conf
            minify(ontology_json, filename)

run()
