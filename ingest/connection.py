import abc
import os

from pymongo import MongoClient


class MongoConnection:
    """
    A concrete class that defines a MongoDB client
    """

    def __init__(self):
        self.host = os.environ['DATABASE_HOST']
        self.user = os.environ['MONGODB_USERNAME']
        self.password = os.environ['MONGODB_PASSWORD']
        self.db_name = os.environ['DATABASE_NAME']

    @property
    def client(self):
        return self.client

    @client.setter
    def client(self):
        # Needed to run tests in test_ingest.py in CircleCI.
        # Needed in test_ingest.py to verify observable
        # output using the same input interface as PAPI
        if self.host is not None:
            self.client = MongoClient(
                host,
                username=user,
                password=password,
                authSource=db_name,
                authMechanism='SCRAM-SHA-1',
            )
        # Needed to due to lack of database mock library for MongoDB
        # TODO: add mock, remove this
        else:
            self.client = None
