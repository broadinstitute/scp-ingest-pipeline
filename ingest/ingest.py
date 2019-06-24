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
From python expression file:
import ingest
ingest.connect(<extract function>, <transform function>)
ingest.ingest()
Ex:
import ingest
ingest_service = ingest.connect(self.extract, self.transform)
ingest_service.ingest()
"""
import argparse
from typing import Union, Generator, Dict, List, Tuple
import time
import os

from google.cloud import firestore
import numpy as np
from gene_data_model import Gene
import pprint
from ingest_loom import Loom
from mtx import Mtx
from dense_matrix import Dense

# Ingest file types

EXPRESSION_FILE_TYPES = ['loom', 'dense', 'mtx']


class IngestService(object):
    def __init__(self, *, matrix_file: str, matrix_file_type: str,
        matrix_bundle: List[str] =None):
        """Initializes variables in ingest service.

        Args:
            extract_fn:
                A function that extracts data for the given file
            transform_fn:
                A function that transforms extracted data into db datamodel

        Returns:
            Nothing
        """
        if not os.path.exists(matrix_file):
		          raise IOError(f"File '{file_path}' not found")
        self.matrix_file_path=  matrix_file
        self.matrix_file_type =  matrix_file_type
        self.matrix_bundle = matrix_bundle
        self.matrix = self.initialize_file_connection()
        self.db = firestore.Client()


    def initialize_file_connection(self):
        file_connections = {
                'loom': Loom,
                'dense': Dense,
                'mtx' : Mtx,
                }

        return(file_connections.get(self.matrix_file_type)(self.matrix_file_path))

    def close_matrix(self):
        self.matrix.close()

    def load_expression_data(self, list_of_transformed_data: List[Gene]) -> None:
        """Loads expression data into firestore

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



    def ingest_expression(self):
        """ Ingests expression files. Calls file type's extract and transform
        functions. Then loads data into firestore.

        Args:
            None

        Returns:
            Nothing
        """
        for data in self.matrix.extract():
            transformed_data = self.matrix.transform_expression_data(*data)
            self.load_expression_data(transformed_data)

def parse_args():
    args = argparse.ArgumentParser(
    prog = 'ingest.py',
    description= __doc__,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    args.add_argument('--matrix-file',
    help='Absolute or relative path to expression file. \
    For 10x data this is the .mtx file',
    )

    subargs = args.add_subparsers()

    ## Ingest Expression files subparser
    parser_ingest_xpression = subargs.add_parser('ingest_expression',
        help='Indicates which expression files are being ingested')

    parser_ingest_xpression.add_argument('--matrix-file-type',
        choices= EXPRESSION_FILE_TYPES, type=str.lower,
        required = True,
        help='Type of expression file that is ingested'
    )
    parser_ingest_xpression.add_argument(
        '--matrix-bundle', default=None, nargs='+',
        help='Names of .genes.tsv, and .barcodes.tsv files'
    )

    parsed_args = args.parse_args()

    if(parsed_args.matrix_file_type == 'mtx'):
        if parsed_args.matrix_bundle is None:
            raise ValueError('Mtx files must include .genes.tsv, and \
            .barcodes.tsv files. See --help for more information')

    return parsed_args

def main():
    """This function handles the actual logic of this script."""

    opts = vars(parse_args())
    print(f'{opts}')
    ingest = IngestService(**opts)

    if hasattr(ingest, 'ingest_expression'):
        getattr(ingest, 'ingest_expression')()



if __name__ == "__main__":
    main()
