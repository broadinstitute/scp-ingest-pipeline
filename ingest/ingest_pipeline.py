"""Ingest Pipeline for expression files and eventually metadata and cluster
files into Firestore.

DESCRIPTION
This cli currently takes in extract and transform functions from different
file types then uploads them into Firestore.

PREREQUISITES
You must have Google Cloud Firestore installed, authenticated, and
configured. Must have Python 3.6 or higher.

EXAMPLES
# Takes expression file and stores it into Firestore

# Ingest dense file
$python ingest_pipeline.py ingest_expression --matrix-file ../tests/data/dense_matrix_19_genes_100k_cells.txt --matrix-file-type dense

# Ingest mtx files
$python ingest_pipeline.py ingest_expression --matrix-file ../tests/data/matrix.mtx --matrix-file-type mtx --gene-file ../tests/data/genes.tsv --barcode-file ../tests/data/barcodes.tsv
"""
import argparse
import os
import time
from typing import Dict, Generator, List, Tuple, Union

import numpy as np
from dense import Dense
from gene_data_model import Gene
from google.api_core import exceptions
from google.cloud import firestore
from mtx import Mtx

# Ingest file types
EXPRESSION_FILE_TYPES = ['dense', 'mtx']


class IngestPipeline(object):
    def __init__(self, *, matrix_file: str, matrix_file_type: str,
                 barcode_file: str = '', gene_file: str = ''):
        """Initializes variables in Ingest Pipeline.

        Args:
            matrix_file: str,
                For expression files, the relative or absolute path to the
                    matrix file
            matrix_file_type: str,
                The matrix file type
            matrix_bundle: List[str]
                Used for MTX files. The matrix bundle consister of the barcode
                    and gene files.

        Returns:
            None
        """
        if matrix_file[:5] != 'gs://' and not os.path.exists(matrix_file):
            raise IOError(f"File '{matrix_file}' not found")
        self.matrix_file_path = matrix_file
        self.matrix_file_type = matrix_file_type
        self.gene_file = gene_file
        self.barcodes_file = barcode_file
        self.matrix = self.initialize_file_connection(matrix_file_type)
        self.db = firestore.Client()

    def initialize_file_connection(self, file_type):
        """Initializes connection to file.

        Args:
            None
        Returns:
            File object.
        """
        # Mtx file types not included because class declaration is different
        file_connections = {
            'dense': Dense
        }
        if self.matrix_file_type == 'mtx':
            return Mtx(self.matrix_file_path, self.gene_file, self.barcodes_file)
        else:
            return(
                file_connections.get(file_type)(self.matrix_file_path))

    def close_matrix(self):
        """Closes connection to file.

        Args:
            None
        Returns:
            None
        """
        self.matrix.close()

    def load_expression_data(self, list_of_expression_models: List[Gene]) -> None:
        """Loads expression data into Firestore.

        Args:
            list_of_transformed_data : List[Gene]
                A list of object type Gene that's stored into Firestore

        Returns:
            None
        """

        # for expression_model in list_of_expression_models:
        for expression_model in list_of_expression_models:
            collection_name = expression_model.get_collection_name()
            doc_ref = self.db.collection(collection_name).document()
            doc_ref.set(expression_model.top_level_doc)
            if expression_model.has_subcollection_data():
                try:
                    subcollection_name = expression_model.get_subcollection_name()
                    doc_ref_sub = doc_ref.collection(
                        subcollection_name).document()
                    doc_ref_sub.set(expression_model.subdocument)
                except exceptions.InvalidArgument as e:
                    # Catches invalid argument exception, which error "Maximum
                    # document size" falls under
                    print(e)
                    batch = self.db.batch()
                    for subdoc in expression_model.chunk_gene_expression_documents():

                        subcollection_name = expression_model.get_subcollection_name()
                        doc_ref_sub = doc_ref.collection(
                            subcollection_name).document()
                        batch.set(doc_ref_sub, subdoc)

                    batch.commit()

    def ingest_expression(self) -> None:
        """Ingests expression files. Calls file type's extract and transform
        functions. Then loads data into Firestore.

        Args:
            None

        Returns:
            None
        """
        if self.gene_file is not None:
            self.matrix.extract()
            transformed_data = self.matrix.transform_expression_data_by_gene()
        else:
            for data in self.matrix.extract():
                transformed_data = self.matrix.transform_expression_data_by_gene(
                    *data)
        self.load_expression_data(transformed_data)
        self.close_matrix()


def create_parser():
    """Creates parser for input arguments.

    Structuring the argument parsing code like this eases automated testing.

    Args:
        None

    Returns:
        parser: ArgumentParser object
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    subparser = parser.add_subparsers()

    # Ingest expression files subparser
    parser_ingest_expression = subparser.add_parser('ingest_expression',
                                                  help='Indicates that expression'
                                                  ' files are being ingested')

    parser_ingest_expression.add_argument('--matrix-file', required=True,
                                          help='Absolute or relative path to '
                                          'expression file. For 10x data this is '
                                          'the .mtx file')

    matrix_file_type_txt = 'Type of expression file that is ingested. If mtx \
        files are being ingested, .genes.tsv and .barcodes.tsv files must be \
        included using --barcode-file <barcode file path> and --gene-file \
        <gene file path>. See --help for more information'

    parser_ingest_expression.add_argument('--matrix-file-type',
                                          choices=EXPRESSION_FILE_TYPES,
                                          type=str.lower,
                                          required=True,
                                          help=matrix_file_type_txt
                                          )

    # Gene and Barcode arguments for MTX bundle
    parser_ingest_expression.add_argument('--barcode-file',
                                          help='Path to .barcodes.tsv files')
    parser_ingest_expression.add_argument('--gene-file',
                                          help='Path to .genes.tsv file')

    parser_ingest_cluster = subparser.add_parser('ingest_clusters',
                                               help='Indicates that cluster'
                                               ' files are being ingested')
    parser_ingest_cluster.add_argument('--cluster-file',
                                       help='Path to cluster files')
    return parser

def validate_arguments(parsed_args):
    """Verify parsed input arguments

    Args:
        parsed_args: Parsed input arguments

    Returns:
        None
    """
    if parsed_args.matrix_file_type == 'mtx' and (parsed_args.gene_file == None
                                                    or parsed_args.barcode_file == None):
        raise ValueError(
            ' Missing arguments: --gene-file and --barcode-file. Mtx files '
            'must include .genes.tsv, and .barcodes.tsv files. See --help for '
            'more information')

def main() -> None:
    """This function handles the actual logic of this script.

    Args:
        None

    Returns:
        None
    """
    
    parsed_args = create_parser().parse_args()
    validate_arguments(parsed_args)
    arguments = vars(parsed_args)
    ingest = IngestPipeline(**arguments)

    if hasattr(ingest, 'ingest_expression'):
        getattr(ingest, 'ingest_expression')()


if __name__ == "__main__":
    main()
