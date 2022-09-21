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
    import random
    import math

    MAX_ATTEMPTS = 5
    # Adopted from https://stackoverflow.com/questions/46939285

    def retry(attempt_num):
        if attempt_num < MAX_ATTEMPTS - 1:
            exp_backoff = pow(2, attempt_num)
            max_jitter = math.ceil(exp_backoff * 0.2)
            final_wait_time = exp_backoff + random.randint(
                0, max_jitter
            )  # exponential back off
            dev_logger.warning(" Waiting %.1f seconds.", final_wait_time)
            time.sleep(final_wait_time)

    @functools.wraps(mongo_op_func)
    def wrapper(*args, **kwargs):
        args = list(args)
        for attempt in range(MAX_ATTEMPTS):
            try:
                return mongo_op_func(*args, **kwargs)
            except AutoReconnect as e:
                if attempt < MAX_ATTEMPTS - 1:
                    dev_logger.warning("PyMongo auto-reconnecting... %s.", str(e))
                    retry(attempt)
                else:
                    raise e
            except BulkWriteError as bwe:
                if attempt < MAX_ATTEMPTS - 1:
                    dev_logger.warning(
                        "Batch ops error occurred. Reinsert attempt %s.", str(attempt)
                    )
                    error = bwe.details["writeErrors"][0]
                    # Check error code to see if any failures are due to violating a unique index (error code 11000)
                    # and discard those documents before retrying
                    if error['code'] == 11000:
                        filtered_docs = discard_inserted_documents(error['op'], args[0], args[1], args[2])
                        if len(filtered_docs) > 0:
                            args[0] = filtered_docs
                            retry(attempt)
                        else:
                            return args[0]
                    else:
                        dev_logger.debug(str(bwe.details))
                        raise bwe
                else:
                    dev_logger.debug(str(bwe.details))
                    raise bwe

    return wrapper


def discard_inserted_documents(error_doc, original_documents, collection_name, mongo_client):
    """Discard any documents that have already been inserted which are violating index constraints
       such documents will have an error code of 11000 for a DuplicateKey error
       from https://github.com/mongodb/mongo/blob/master/src/mongo/base/error_codes.yml#L467
       Uses Mongo query to get existing names and filters against documents to be inserted as
       only first error is returned, not all possible errors

       Parameters:
           error_doc (Dict): document that failed to insert in original transaction
           original_documents (List[Dict]): list of documents from original transaction that failed
           collection_name (String): name of collection where documents were to be inserted
           mongo_client (MongoClient): MongoDB client to perform query for existing documents

       Returns:
           List[Dict]: list of documents with existing entries removed
    """
    names_to_find = list(doc['name'] for doc in original_documents)
    query = {
        'study_id' : error_doc['study_id'],
        'study_file_id' : error_doc['study_file_id'],
        'name' : { "$in" : names_to_find }
    }
    if 'linear_data_type' in error_doc:
        query['linear_data_type'] = error_doc['linear_data_type']
    fields = { 'name': 1 }
    raw_names = mongo_client[collection_name].find(query, fields)
    existing_names = list(doc['name'] for doc in raw_names)
    return list(doc for doc in original_documents if doc['name'] not in existing_names)
