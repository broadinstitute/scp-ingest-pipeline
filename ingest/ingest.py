from typing import Union, Generator, Dict, List, Tuple
import time

from google.cloud import firestore
import numpy as np
from gene_data_model import Gene
import pprint

class IngestService:
    def __init__(self,extract_fn, transform_fn):
        """Initializes variables in ingest service.

        Args:
            extract_fn:  A function that extracts data for the given file
            transform_fn: A function that transforms extracted data into db
                datamodel

        Returns:
            Nothing
        """
        self.extract = extract_fn
        self.transform = transform_fn
        self.db = firestore.Client()


    def load(self, list_of_transformed_data: List[Gene]) -> None:
        """Loads data into firestore

        Args:
            list_of_transformed_data (List[Gene]): A list of object type Gene
                that's stored into Firestore

        Returns:
            Nothing
        """
        batch = self.db.batch()
        for transformed_data in list_of_transformed_data:
            for collection, document in transformed_data.items():
                for gene, data in document.items():
                    gene_doc_ref = self.db.collection(collection).document(gene)
                    batch.set(gene_doc_ref, data)
                    batch.commit()
                    time.sleep(.2)


    def ingest(self):
        """
        """
        for data in self.extract():
            transformed_data = self.transform(*data)
            self.load(transformed_data)


def connect(extract_fn, transform_fn) -> IngestService:
    """
    Connects to Ingest service.

    Args:
        extract_fn :  A function that extracts data for the given file
        transform_fn: A function that transforms extracted data into db datamodel


    Returns:
        An Ingest_Service instance
    """
    return IngestService(extract_fn, transform_fn)
