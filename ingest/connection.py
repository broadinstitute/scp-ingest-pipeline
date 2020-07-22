import abc
import os

from pymongo import MongoClient


class MongoConnection:
    """
    A concrete class that defines a MongoDB client
    """

    def __init__(self):
        # Needed to run tests in test_ingest.py in CircleCI.
        # Needed in test_ingest.py to verify observable
        # output using the same input interface as PAPI
        if os.environ.get('DATABASE_HOST') is not None:
            self.host = os.environ['DATABASE_HOST']
            self.user = os.environ['MONGODB_USERNAME']
            self.password = os.environ['MONGODB_PASSWORD']
            self.db_name = os.environ['DATABASE_NAME']
            client = MongoClient(
                self.host,
                username=self.user,
                password=self.password,
                authSource=self.db_name,
                authMechanism='SCRAM-SHA-1',
            )
            self._client=client[self.db_name]
        # Needed to due to lack of database mock library for MongoDB
        # TODO: add mock, remove this
        else:
            self._client = None
