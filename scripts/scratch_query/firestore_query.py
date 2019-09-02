"""Firestore query examples using python client.

DESCRIPTION
This cli runs a fixed set of queries against the Firestore cell_metadata collection.

PREREQUISITES
You must have Google Cloud Firestore installed, authenticated, and
configured. Must have Python 3.6 or higher. Indexing must be turned off for
all collections.

EXAMPLES
# Query Firestore with fixed queries
python firestore_query.py is_living array_contains yes
"""
import argparse
from google.cloud import firestore
import warnings

warnings.filterwarnings(
    "ignore", "Your application has authenticated using end user credentials"
)


def create_parser():
    """Parse command line values for firestore_query
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
    parser.add_argument('metadata', help='study metadata being queried')
    parser.add_argument('comparison_operation', help='comparison operation for query')
    parser.add_argument('value', help='value for comparison')
    return parser


def match_query(query_params):
    metadata, operator, value = query_params
    matching_studies = set()
    db = firestore.Client()
    cm = db.collection(u'cell_metadata')
    living_query = cm.where(u'name', u'==', metadata).stream()
    for doc in living_query:
        query_path = 'cell_metadata/' + doc.id + '/data'
        # foo = "cell_metadata/jkIaZbCKXdyEeZOhSDGS/data"
        # subdocs = db.collection(query_path).stream()  # all subdocs
        subdocs = db.collection(query_path).where(u'values', operator, value).stream()
        # print(u'{} => {}'.format(doc.study_accession))
        # results = subdocs.stream()
        for subdoc in subdocs:
            # print(u'{} => {}'.format(subdoc.id, subdoc.to_dict()))
            test = doc.to_dict()
            matching_studies.add(test['study_accession'])
            # for key, value in doc.to_dict().items():
            #     print(key, value)
            # print(u'{} => {}'.format(doc.study_accession))
    if matching_studies:
        return matching_studies
    else:
        return None


if __name__ == '__main__':
    args = create_parser().parse_args()
    query_params = (args.metadata, args.comparison_operation, args.value)
    matches = match_query(query_params)
    print('Studies matching query:', matches)
