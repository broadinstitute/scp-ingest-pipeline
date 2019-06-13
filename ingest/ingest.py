from typing import Union, Generator, Dict, List, Tuple

from google.cloud import firestore
import numpy as np
from gene_data_model import Gene
import pprint

class IngestService:
    def __init__(self, file_path, file_type_class):
        """Initializes variables in ingest service.

        Args:
            extract_fn:  A function that extracts data for the given file
            transform_fn: A function that transforms extracted data into db
                datamodel

        Returns:
            Nothing
        """
        if not os.path.exists(file_path):
		          raise IOError(f"File '{file_path}' not found")
        file_type_class = file_type_class(file_path)
        self.extract = file_type_class.extraxt
        self.transform = file_type_class.transform
        self.db = firestore.Client()


    def load(self, list_of_transformed_data: List[Gene]) -> None:
        """
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog ='ingest.py',
        description= __doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    #Positional argument
    parser.add_argument(
        "file-path",
        help='Absolute or relative path to expression file'
    )


    subargs = args.add_subparsers()

    parser_ingest_expression_file = subargs.add_parser('ingest-expression-files',
        dest='ingest_expression_files',
        description = "Indicates that expression files will be ingested")

    parser_ingest_expression_file.add_argument(
        '--file-type',
        required=True,
        choices=['loom', 'dense-matrix', 'mtx'],
        help = 'Indicates what type of expression file willbe ingest'
    )

    dense_mtx = subargs.add_parser('dense_matrx', help='Ingest expression file that is a dense matrix')
    dense_mtx.add_argument('--dense-matrix', help = 'Indicates that expression file is a dense matrix file')

    mtx = subargs.add_parser('mtx', help='Ingest expression file that is a mtx')
    mtx.add_argument('mtx')
    args = parser.parse_args()

    parsed_args = args.parse_args()

    if hasarrt(parsed_args.file_type, 'loom') :
        ingest_loom = Ingest_Service(extract_fn, transform_fn)

def connect(extract_fn, transform_fn) -> Ingest_Service:
    """
    Connects to Ingest service.

    Args:
        extract_fn :  A function that extracts data for the given file
        transform_fn: A function that transforms extracted data into db datamodel


    Returns:
        An Ingest_Service instance
    """
    return Ingest_Service(extract_fn, transform_fn)
