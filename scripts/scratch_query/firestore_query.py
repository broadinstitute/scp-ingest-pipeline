"""Firestore query examples using python client.

DESCRIPTION
This cli runs a fixed set of queries against the Firestore cell_metadata collection.

PREREQUISITES
You must have Google Cloud Firestore installed, authenticated, and
configured. Must have Python 3.6 or higher. Indexing must be turned off for
all collections.

EXAMPLES
# Query Firestore with fixed queries
python firestore_query.py
"""

from google.cloud import firestore

matching_studies = set()
db = firestore.Client()
cm = db.collection(u'cell_metadata')
living_query = cm.where(u'name', u'==', u'is_living').stream()
for doc in living_query:
    query_path = 'cell_metadata/' + doc.id + '/data'
    # foo = "cell_metadata/jkIaZbCKXdyEeZOhSDGS/data"
    # subdocs = db.collection(query_path).stream()  # all subdocs
    subdocs = (
        db.collection(query_path).where(u'values', u'array_contains', u'yes').stream()
    )
    # print(u'{} => {}'.format(doc.study_accession))
    # results = subdocs.stream()
    for subdoc in subdocs:
        # print(u'{} => {}'.format(subdoc.id, subdoc.to_dict()))
        test = doc.to_dict()
        matching_studies.add(test['study_accession'])
        # for key, value in doc.to_dict().items():
        #     print(key, value)
        # print(u'{} => {}'.format(doc.study_accession))
print('Studies matching query:', matching_studies)

# yes_living_query = (
#     db.collection(u'cell_metadata/jkIaZbCKXdyEeZOhSDGS/data').where(u'values', u'==', u'yes').stream()
# )
# for doc in yes_living_query:
#     print(u'{} => {}'.format(doc.id, doc.to_dict()))


# cm = db.collection(u'cell_metadata')
# all_docs = cm.stream()
# for doc in all_docs:
#     print(u'{} => {}'.format(doc.id, doc.to_dict()))

# living_query = cm.where(u'name', u'==', u'is_living')

# for doc in living_query:
#     print(u'{} => {}'.format(doc.id, doc.to_dict()))
