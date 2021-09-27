import unittest
import sys
sys.path.append("../ingest")
from mongo_connection import discard_inserted_documents


class TestDiscardDocuments(unittest.TestCase):
    # test discarding documents from insert_many callback that have 11000 error code
    def test_discard_inserted_documents(self):
        error_docs = [
            {
                'code': 11000,
                'op': {'id': 1, 'foo': 'bar', 'bing': 'baz'}
            },
            {
                'code': 11000,
                'op': {'id': 2, 'foo': 'bar', 'bing': 'blurg'}
            },
            {
                'code': 11001,
                'op': {'id': 3, 'blah': 'blurg', 'boo': 'bah'}
            }
        ]

        original_docs = [
            {'id': 1, 'foo': 'bar', 'bing': 'baz'},
            {'id': 2, 'foo': 'bar', 'bing': 'blurg'},
            {'id': 3, 'blah': 'blurg', 'boo': 'bah'},
            {'id': 4, 'blah': 'barf', 'boo': 'bee'}
        ]

        filtered_docs = discard_inserted_documents(error_docs, original_docs)
        self.assertEqual(2, len(filtered_docs), 'Did not correctly return 2 filtered documents')
        ids = list(doc['id'] for doc in filtered_docs)
        ids.sort()
        self.assertEqual([3, 4], ids, 'Did not return correct IDs')
