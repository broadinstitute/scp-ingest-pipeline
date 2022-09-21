import unittest
import sys
sys.path.append("../ingest")
from mongo_connection import discard_inserted_documents
from unittest.mock import MagicMock

client_values = { "data_arrays": MagicMock() }
client_values["data_arrays"].find.return_value = [
    { 'id': 1, 'name': 'foo' },
    { 'id': 2, 'name': 'bar' }
]

class TestDiscardDocuments(unittest.TestCase):
    # mock MongoDB client to return data from client_values
    client_mock: MagicMock = MagicMock()
    client_mock.keys.return_value.__iter__.return_value = client_values.keys()
    client_mock.__getitem__.side_effect = lambda k: client_values[k]

    # test discarding documents from insert_many callback that have 11000 error code
    def test_discard_inserted_documents(self):
        error_doc = { 'id': 1, 'name': 'foo', 'study_id': 1, 'study_file_id': 1, 'linear_data_type': 'Gene' }
        original_docs = [
            {'id': 1, 'name': 'foo', 'study_id': 1, 'study_file_id': 1, 'linear_data_type': 'Gene'},
            {'id': 2, 'name': 'bar', 'study_id': 1, 'study_file_id': 1, 'linear_data_type': 'Gene'},
            {'id': 3, 'name': 'blurg', 'study_id': 1, 'study_file_id': 1, 'linear_data_type': 'Gene'},
            {'id': 4, 'name': 'barf', 'study_id': 1, 'study_file_id': 1, 'linear_data_type': 'Gene'}
        ]

        filtered_docs = discard_inserted_documents(error_doc,
                                                   original_docs,
                                                   'data_arrays',
                                                   TestDiscardDocuments.client_mock)
        self.assertEqual(2, len(filtered_docs), 'Did not correctly return 2 filtered documents')
        ids = list(doc['id'] for doc in filtered_docs)
        ids.sort()
        self.assertEqual([3, 4], ids, 'Did not return correct IDs')
