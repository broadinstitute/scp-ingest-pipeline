import os
import functools
import logging
import time
from pymongo import MongoClient
from pymongo.errors import AutoReconnect, BulkWriteError


class MongoConnection:
    """
    A concrete class that defines a MongoDB client
    """

    MAX_AUTO_RECONNECT_ATTEMPTS = 5

    def __init__(self):
        # Needed to run tests in test_ingest.py in CircleCI.
        # Needed in test_ingest.py to verify observable
        # output using the same input interface as PAPI
        if os.environ.get("DATABASE_HOST") is not None:
            self.host = os.environ["DATABASE_HOST"]
            self.user = os.environ["MONGODB_USERNAME"]
            self.password = os.environ["MONGODB_PASSWORD"]
            self.db_name = os.environ["DATABASE_NAME"]
            client = MongoClient(
                self.host,
                username=self.user,
                password=self.password,
                authSource=self.db_name,
                authMechanism="SCRAM-SHA-1",
            )
            self._client = client[self.db_name]
        # Needed to due to lack of database mock library for MongoDB
        # TODO: add mock, remove this
        else:
            self._client = None


# Adopted from https://gist.github.com/anthonywu/1696591#file
# -graceful_auto_reconnect-py
def graceful_auto_reconnect(mongo_op_func):
    """Gracefully handles a reconnection event as well as other exceptions
        for mongo.
    """
    MAX_ATTEMPTS = 5

    def retry(attempt_num):
        if attempt_num < MAX_ATTEMPTS - 1:
            wait_time = 0.5 * pow(2, attempt_num)  # exponential back off
            logging.warning(" Waiting %.1f seconds.", wait_time)
            time.sleep(wait_time)

    @functools.wraps(mongo_op_func)
    def wrapper(*args, **kwargs):
        for attempt in range(MAX_ATTEMPTS):
            try:
                return mongo_op_func(*args, **kwargs)
            except AutoReconnect as e:
                if attempt < MAX_ATTEMPTS - 1:
                    logging.warning("PyMongo auto-reconnecting... %s.", str(e))
                    retry(attempt)
                else:
                    raise e
            except BulkWriteError as bwe:
                if attempt < MAX_ATTEMPTS - 1:
                    logging.warning(
                        "Batch ops error occurred. Reinsert attempt %s.", str(attempt)
                    )
                    retry(attempt)
                else:
                    raise BulkWriteError(bwe.details)
            except Exception as e:
                raise e

    return wrapper
