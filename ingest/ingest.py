"""
Ingest Service for expression files and extentually metadata and cluster files
into firestore.

DESCRIPTION
This module currently takes in extract and transform functions from file types
then uploads them into Firestore.

PREREQUISITES
You must have google could Firestore installed, authenticated
 configured. Must have python 3.6 or higher.

EXAMPLES
# Takes expression file and stores it into firestore
From expression file:
import ingest
ingest.connect(<extract function>, <transform function>)
ingest.ingest()
Ex:
import ingest
ingest_service = ingest.connect(self.extract, self.transform)
ingest_service.ingest()
"""

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
            extract_fn:
                A function that extracts data for the given file
            transform_fn:
                A function that transforms extracted data into db datamodel

        Returns:
            Nothing
        """
        self.extract = extract_fn
        self.transform = transform_fn
        self.db = firestore.Client()


    def load(self, list_of_transformed_data: List[Gene]) -> None:
        """Loads data into firestore

        Args:
            list_of_transformed_data : List[Gene]
                A list of object type Gene that's stored into Firestore

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
        """ Ingests files. Calls provided extract and transform functions. Then
         loads data into firestore.

        Args:
            None

        Returns:
            Nothing
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
