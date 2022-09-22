import unittest
import sys
sys.path.append("../ingest")
from mongo_connection import discard_inserted_documents

class TestDiscardDocuments(unittest.TestCase):
    # test discarding documents from insert_many callback that have 11000 error code
    def test_discard_inserted_documents(self):
        error_docs = [
            {'code': 11000, 'op': {'id': 1, 'name': 'foo', 'array_index': 0}},
            {'code': 11000, 'op': {'id': 3, 'name': 'bar', 'array_index': 0}},
        ]
        original_docs = [
            {'id': 1, 'name': 'foo', 'array_index': 0},
            {'id': 2, 'name': 'foo', 'array_index': 1},
            {'id': 3, 'name': 'bar', 'array_index': 0},
            {'id': 4, 'name': 'blurg','array_index': 0},
            {'id': 5, 'name': 'barf', 'array_index': 0}
        ]
        filtered_docs = discard_inserted_documents(error_docs, original_docs)
        self.assertEqual(3, len(filtered_docs), 'Did not correctly return 3 filtered documents')
        ids = list(doc['id'] for doc in filtered_docs)
        ids.sort()
        self.assertEqual([2, 4, 5], ids, 'Did not return correct IDs')
