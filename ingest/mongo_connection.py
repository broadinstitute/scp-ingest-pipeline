from __future__ import annotations
import os
import functools
import time
from pymongo import MongoClient
from pymongo.errors import AutoReconnect, BulkWriteError

try:
    from monitor import setup_logger
except ImportError:
    from .monitor import setup_logger
dev_logger = setup_logger(__name__, "log.txt", format="support_configs")


class MongoConnection:
    """
    A concrete class that defines a MongoDB client
    """

    MAX_AUTO_RECONNECT_ATTEMPTS = 10

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
                serverSelectionTimeoutMS=60_000
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
    import random
    import math

    # Adopted from https://stackoverflow.com/questions/46939285

    def retry(attempt_num):
        if attempt_num < MongoConnection.MAX_AUTO_RECONNECT_ATTEMPTS - 1:
            exp_backoff = pow(2, attempt_num)
            max_jitter = math.ceil(exp_backoff * 0.5)
            final_wait_time = exp_backoff + random.randint(
                0, max_jitter
            )  # exponential back off
            dev_logger.warning(" Waiting %.1f seconds.", final_wait_time)
            time.sleep(final_wait_time)

    @functools.wraps(mongo_op_func)
    def wrapper(*args, **kwargs):
        args = list(args)
        for attempt in range(MongoConnection.MAX_AUTO_RECONNECT_ATTEMPTS):
            try:
                return mongo_op_func(*args, **kwargs)
            except AutoReconnect as e:
                if attempt < MongoConnection.MAX_AUTO_RECONNECT_ATTEMPTS - 1:
                    dev_logger.warning("PyMongo auto-reconnecting... %s.", str(e))
                    retry(attempt)
                else:
                    raise e
            except BulkWriteError as bwe:
                if attempt < MongoConnection.MAX_AUTO_RECONNECT_ATTEMPTS - 1:
                    dev_logger.warning(
                        "Batch ops error occurred. Reinsert attempt %s.", str(attempt)
                    )
                    error = bwe.details["writeErrors"]
                    # Check error code to see if any failures are due to violating a unique index (error code 11000)
                    # and discard those documents before retrying
                    filtered_docs = discard_inserted_documents(error, args[0])
                    if len(filtered_docs) > 0:
                        args[0] = filtered_docs
                        retry(attempt)
                    else:
                        return args[0]
                else:
                    log_error_without_values(bwe)
                    raise bwe

    return wrapper


def log_error_without_values(error):
    """Remove 'values' array from log messages as this is usually very long and makes logs difficult to read

    :param error: BulkWriteError from failed transaction
    """
    for entry in error.details["writeErrors"]:
        if 'values' in entry['op']:
            entry['op']['values'] = '<filtered>'
    dev_logger.debug(str(error.details))


def discard_inserted_documents(errors, original_documents) -> list[dict]:
    """Discard any documents that have already been inserted which are violating index constraints
    such documents will have an error code of 11000 for a DuplicateKey error
    from https://github.com/mongodb/mongo/blob/master/src/mongo/base/error_codes.yml#L467

    :param errors: (list[dict]) Documents that failed to insert from transaction
    :param original_documents: (list[dict]) list of documents from original transaction that failed
    :returns list[dict]: list of documents with existing entries removed
    """
    error_docs = []
    for doc in errors:
        if doc['code'] == 11000:
            error_docs.append(
                {'name': doc['op']['name'], 'array_index': doc['op']['array_index']}
            )
    return list(
        d
        for d in original_documents
        if {'name': d['name'], 'array_index': d['array_index']} not in error_docs
    )
