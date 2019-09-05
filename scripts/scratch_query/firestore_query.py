"""Firestore query examples using python client.

DESCRIPTION
This cli runs a fixed set of queries against the Firestore cell_metadata collection.

PREREQUISITES
You must have Google Cloud Firestore installed, authenticated, and
configured. Must have Python 3.6 or higher. Indexing must be turned off for
all collections.

EXAMPLES
# Query Firestore cell_metadata, supplying known good metadata names and values queries
python firestore_query.py --query is_living,yes
python firestore_query.py --query species__ontology_label,"Homo sapiens"
python firestore_query.py --query species__ontology_label,"Homo sapiens" --query is_living,yes
python firestore_query.py --query species__ontology_label,"Homo sapiens" --query sex,female
python firestore_query.py --query species__ontology_label,"Homo sapiens" --query sex,male
python firestore_query.py --query species__ontology_label,"Homo sapiens" --query sex,unknown
python firestore_query.py --query species__ontology_label,"Homo sapiens" --query is_living,yes --query sex,female
python firestore_query.py --query diagnosis,FGID
"""
import argparse
from google.cloud import firestore
import warnings
import sys

warnings.filterwarnings(
    "ignore", "Your application has authenticated using end user credentials"
)


def create_parser():
    """Parse command line values for firestore_query
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    def to_list(arg):
        return [i for i in arg.split(",")]

    parser.add_argument(
        '--query', help='metadata,value being queried', action="append", type=to_list
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser


def match_query(query_params):
    metadata, value = query_params
    matching_studies = set()
    db = firestore.Client()
    cm = db.collection(u'cell_metadata')
    living_query = cm.where(u'name', u'==', metadata).stream()
    for doc in living_query:
        query_path = 'cell_metadata/' + doc.id + '/data'
        subdocs = (
            db.collection(query_path)
            .where(u'values', u'array_contains', value)
            .stream()
        )
        for subdoc in subdocs:
            # print(u'{} => {}'.format(subdoc.id, subdoc.to_dict()))
            test = doc.to_dict()
            matching_studies.add(test['study_accession'])
    if matching_studies:
        return matching_studies
    else:
        return None


if __name__ == '__main__':
    args = create_parser().parse_args()
    results = []
    for query in args.query:
        matches = match_query(query)
        metadata, value = query
        print('Studies with cells where', metadata, 'is', value, ':', matches)
        results.append(matches)
    if len(args.query) > 1:
        print('++++++++++++++++++++')
        try:
            print(
                'Studies with cells matching all criteria:', set.intersection(*results)
            )
        except TypeError:
            print('No matching studies - an input query returned no matches')
