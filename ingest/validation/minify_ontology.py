import json

def get_synonyms(node):
    if 'meta' not in node or 'synonyms' not in node['meta']:
        return ''

    raw_synonyms = []
    synonym_nodes = node['meta']['synonyms']
    for synonym_node in synonym_nodes:
        raw_synonyms.append(synonym_node['val'])
    print('raw_synonyms', raw_synonyms)
    synonyms = '||'.join(raw_synonyms) # Unambiguously delimit synonyms
    return synonyms


def minify(ontology_shortname):
    ontology_shortname_lc = ontology_shortname.lower()
    with open(f'{ontology_shortname_lc}.json') as f:
        ontology = json.loads(f.read())

    nodes = ontology['graphs'][0]['nodes']

    raw_disease_nodes = list(filter(lambda n: f'obo/{ontology_shortname}' in n['id'] and 'lbl' in n, nodes))
    disease_nodes = list(map(
        lambda n: '\t'.join(
            [n['id'].split('/')[-1], n['lbl'], get_synonyms(n)]
        ), raw_disease_nodes
    ))

    with open(f'{ontology_shortname_lc}.min.tsv', 'w') as f:
        f.write('\n'.join(disease_nodes))

minify('MONDO')
