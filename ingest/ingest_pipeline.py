"""Ingest Pipeline for ingesting expression, metadata and cluster
files into Firestore.

DESCRIPTION
This cli currently takes in extract and transform functions from different
file types then uploads them into Firestore.

PREREQUISITES
You must have Google Cloud Firestore installed, authenticated, and
configured. Must have Python 3.6 or higher. Indexing must be turned off for sub-collections.

EXAMPLES
# Takes expression file and stores it into Firestore

# Ingest cluster file
python ingest_pipeline.py --study-accession SCP1 --file-id 123abc ingest_cluster --cluster-file ../tests/data/10k_cells_29k_genes.cluster.txt --ingest-cluster --name cluster1 --domain-ranges "{'x':[-1, 1], 'y':[-1, 1], 'z':[-1, 1]}"

# Ingest Cell Metadata file
python ingest_pipeline.py --study-accession SCP1 --file-id 123abc ingest_cell_metadata --cell-metadata-file ../tests/data/10k_cells_29k_genes.metadata.tsv --ingest-cell-metadata

# Ingest dense file
python ingest_pipeline.py  --study-accession SCP1 --file-id 123abc ingest_expression --taxon-name 'Homo sapiens' --taxon-common-name human --ncbi-taxid 9606 --matrix-file ../tests/data/dense_matrix_19_genes_100k_cells.txt --matrix-file-type dense

# Ingest loom file
python ingest_pipeline.py  --study-accession SCP1 --file-id 123abc ingest_expression --matrix-file ../tests/data/test_loom.loom  --matrix-file-type loom --taxon-name 'Homo Sapiens' --taxon-common-name humans

# Subsample cluster and metadata file
python ingest_pipeline.py --study-accession SCP1 --file-id 123abc ingest_subsample --cluster-file ../tests/data/test_1k_cluster_Data.csv --cell-metadata-file ../tests/data/test_1k_metadata_Data.csv --subsample

# Ingest mtx files
python ingest_pipeline.py --study-accession SCP1 --file-id 123abc ingest_expression --taxon-name 'Homo Sapiens' --taxon-common-name humans --matrix-file ../tests/data/matrix.mtx --matrix-file-type mtx --gene-file ../tests/data/genes.tsv --barcode-file ../tests/data/barcodes.tsv
"""
import argparse
from typing import Dict, Generator, List, Tuple, Union  # noqa: F401
import ast
import os
import sys

from cell_metadata import CellMetadata
from clusters import Clusters
from dense import Dense
from gene_data_model import Gene
from google.api_core import exceptions
from google.cloud import firestore
from mtx import Mtx
from subsample import SubSample
from loom import Loom

# Ingest file types
EXPRESSION_FILE_TYPES = ["dense", "mtx", "loom"]


class IngestPipeline(object):
    def __init__(
        self,
        *,
        file_id: str,
        study_accession: str,
        matrix_file: str = None,
        matrix_file_type: str = None,
        cell_metadata_file: str = None,
        cluster_file: str = None,
        subsample=False,
        ingest_cell_metadata=False,
        ingest_cluster=False,
        **kwargs,
    ):
        """Initializes variables in ingest service."""
        self.file_id = file_id
        self.study_accession = study_accession
        self.matrix_file = matrix_file
        self.matrix_file_type = matrix_file_type
        self.db = firestore.Client()
        self.cluster_file = cluster_file
        self.kwargs = kwargs
        self.cell_metadata_file = cell_metadata_file
        if matrix_file is not None:
            self.matrix = self.initialize_file_connection(matrix_file_type, matrix_file)
        elif ingest_cell_metadata:
            self.cell_metadata = self.initialize_file_connection(
                "cell_metadata", cell_metadata_file
            )
        elif ingest_cluster:
            self.cluster = self.initialize_file_connection("cluster", cluster_file)
        elif matrix_file is None:
            self.matrix = matrix_file

    def initialize_file_connection(self, file_type, file_path):
        """Initializes connection to file.

            Returns:
                File object.
        """
        # Mtx file types not included because class declaration is different
        file_connections = {
            "dense": Dense,
            "cell_metadata": CellMetadata,
            "cluster": Clusters,
            "mtx": Mtx,
            "loom": Loom,
        }
        return file_connections.get(file_type)(
            file_path, self.file_id, self.study_accession, **self.kwargs
        )

    def close_matrix(self):
        """Closes connection to file"""
        self.matrix.close()

    def load_expression_data(self, list_of_expression_models: List[Gene]) -> None:
        """Loads expression data into Firestore.

    Args:
        list_of_transformed_data: List[Gene]
           A list of object type Gene that's stored into Firestore

    Returns:
        None
    """
        for expression_model in list_of_expression_models:
            collection_name = expression_model.COLLECTION_NAME
            batch = self.db.batch()
            doc_ref = self.db.collection(collection_name).document()
            batch.set(doc_ref, expression_model.top_level_doc)
            batch.commit()
            i = 0
            if expression_model.has_subcollection_data():
                try:
                    print(f"Ingesting {expression_model.name}")
                    subcollection_name = expression_model.SUBCOLLECTION_NAME
                    doc_ref_sub = doc_ref.collection(subcollection_name).document()
                    print(
                        f"Length of scores is: {len(expression_model.expression_scores)}"
                    )
                    doc_ref_sub.set(expression_model.subdocument)
                except exceptions.InvalidArgument as e:
                    # Catches invalid argument exception, which error "Maximum
                    # document size" falls under
                    print(e)
                    batch = self.db.batch()
                    for subdoc in expression_model.chunk_gene_expression_documents(
                        doc_ref_sub.id, doc_ref_sub._document_path
                    ):
                        print({i})
                        batch.set(doc_ref_sub, subdoc)
                        i += 1

                    batch.commit()

    def load_cell_metadata(self):
        """Loads cell metadata files into firestore."""

        collection_name = self.cell_metadata.COLLECTION_NAME
        subcollection_name = self.cell_metadata.SUBCOLLECTION_NAME
        for annotation in self.cell_metadata.top_level_doc.keys():
            doc_ref = self.db.collection(collection_name).document()
            self.cell_metadata.update_unqiue_values(annotation)
            doc_ref.set(self.cell_metadata.top_level_doc[annotation])
            try:
                print(f"Ingesting {annotation}")
                print(
                    f"Length of values are: {len(self.cell_metadata.data_subcollection[annotation]['values'])}"
                )
                subcollection_doc = self.cell_metadata.data_subcollection[annotation]
                doc_ref_sub = doc_ref.collection(subcollection_name).document()
                doc_ref_sub.set(subcollection_doc)
            except exceptions.InvalidArgument:
                # Catches invalid argument exception, which error "Maximum
                # document size" falls under
                batch = self.db.batch()
                for subdoc in self.cell_metadata.chunk_subdocuments(
                    doc_ref_sub.id, doc_ref_sub._document_path, annotation
                ):
                    batch.set(doc_ref_sub, subdoc)

                batch.commit()
            else:
                return 1
        return 0

    def load_cluster_files(self):
        """Loads cluster files into Firestore."""
        collection_name = self.cluster.COLLECTION_NAME
        doc_ref = self.db.collection(collection_name).document()
        doc_ref.set(self.cluster.top_level_doc)
        subcollection_name = self.cluster.SUBCOLLECTION_NAME
        for annot_name in self.cluster.cluster_subdocs.keys():
            # batch = self.db.batch()
            doc_ref_sub = doc_ref.collection(subcollection_name).document()
            doc_ref_sub.set(self.cluster.cluster_subdocs[annot_name])

    def ingest_expression(self) -> None:
        """Ingests expression files. Calls file type's extract and transform
    functions. Then loads data into Firestore.

    Args:
        None

    Returns:
        None
    """
        transformed_data = []
        if self.kwargs["gene_file"] is not None:
            self.matrix.extract()
            transformed_data = self.matrix.transform_expression_data_by_gene()
        elif self.matrix_file_type == "loom":
            for expression_ds in self.matrix.extract():
                transformed_data = (
                    transformed_data
                    + self.matrix.transform_expression_data_by_gene(expression_ds)
                )
        else:
            while True:
                row = self.matrix.extract()
                if row is None:
                    break
                transformed_data.append(
                    self.matrix.transform_expression_data_by_gene(row)
                )
        self.load_expression_data(transformed_data)
        # # self.close_matrix()

    def ingest_cell_metadata(self):
        """Ingests cell metadata files into Firestore."""
        if self.cell_metadata.is_valid_file:
            print("valid file")
            while True:
                row = self.cell_metadata.extract()
                if row is None:
                    break
                self.cell_metadata.transform(row)
            load_status = self.load_cell_metadata()
            return load_status
        else:
            print("invalid file")
            return 1

    def ingest_cluster(self):
        """Ingests cluster files into Firestore."""
        while True:
            row = self.cluster.extract()
            if row is None:
                self.cluster.update_points()
                self.cluster.update_cell_annotations_field()
                break
            self.cluster.transform(row)
        self.load_cluster_files()

    def subsample(self):
        """Method for subsampling cluster and metadata files"""

        subsample = SubSample(
            cluster_file=self.cluster_file, cell_metadata_file=self.cell_metadata_file
        )

        def create_cluster_subdoc(scope):
            for subdoc in subsample.subsample():
                print(subdoc)
                annot_name = subdoc[1][0]
                annot_type = subdoc[1][1]
                sample_size = subdoc[2]
                for key_value in subdoc[0].items():
                    Clusters.create_cluster_subdoc(
                        key_value[0],
                        annot_type,
                        value=key_value[1],
                        subsample_annotation=f"{annot_name}--{annot_type}--{scope}",
                        subsample_threshold=sample_size,
                    )

        create_cluster_subdoc("cluster")
        if self.cell_metadata_file is not None:
            subsample.prepare_cell_metadata()
            subsample.determine_coordinates_and_cell_names()
            create_cluster_subdoc("study")


def create_parser():
    """Creates parser for input arguments.

    Structuring the argument parsing code like this eases automated testing.

    Args:
        None

    Returns:
        parser: ArgumentParser object
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--study-accession",
        required=True,
        help="Single study accession associated with ingest files",
    )
    parser.add_argument("--file-id", required=True, help="MongoDB identifier")

    subparsers = parser.add_subparsers()

    # Ingest expression files subparsers
    parser_ingest_expression = subparsers.add_parser(
        "ingest_expression",
        help="Indicates that expression" " files are being ingested",
    )

    parser_ingest_expression.add_argument(
        "--matrix-file",
        required=True,
        help="Absolute or relative path to "
        "expression file. For 10x data this is "
        "the .mtx file",
    )

    matrix_file_type_txt = "Type of expression file that is ingested. If mtx \
        files are being ingested, .genes.tsv and .barcodes.tsv files must be \
        included using --barcode-file <barcode file path> and --gene-file \
        <gene file path>. See --help for more information"

    parser_ingest_expression.add_argument(
        "--taxon-name",
        help="Scientific name of taxon associated with file.  E.g. 'Homo sapiens'",
    )
    parser_ingest_expression.add_argument(
        "--taxon-common-name",
        help="Common name of taxon associated with file.  E.g. 'human'",
    )
    parser_ingest_expression.add_argument(
        "--ncbi-taxid",
        help="NCBI Taxonomy ID of taxon associated with file.  E.g. 9606",
    )
    parser_ingest_expression.add_argument(
        "--genome-assembly-accession",
        help="Genome assembly accession for file.  E.g. 'GCA_000001405.15'",
    )
    parser_ingest_expression.add_argument(
        "--genome-annotation",
        help="Genomic annotation for expression files.  E.g. 'Ensembl 94'",
    )

    parser_ingest_expression.add_argument(
        "--matrix-file-type",
        choices=EXPRESSION_FILE_TYPES,
        type=str.lower,
        required=True,
        help=matrix_file_type_txt,
    )

    # Gene and Barcode arguments for MTX bundle
    parser_ingest_expression.add_argument(
        "--barcode-file", help="Path to .barcodes.tsv files"
    )
    parser_ingest_expression.add_argument("--gene-file", help="Path to .genes.tsv file")

    # Parser ingesting cell metadata files
    parser_cell_metadata = subparsers.add_parser(
        "ingest_cell_metadata",
        help="Indicates that cell metadata files are being " "ingested",
    )
    parser_cell_metadata.add_argument(
        "--cell-metadata-file",
        required=True,
        help="Absolute or relative path to cell metadata file.",
    )
    parser_cell_metadata.add_argument(
        "--ingest-cell-metadata",
        required=True,
        action="store_true",
        help="Indicates that subsampliing functionality should be invoked",
    )

    parser_cell_metadata.add_argument(
        "--validate-cell-metadata",
        "-vcm",
        action="store_true",
        help="Indicates that file should be validated regardless on ingest being invoked",
    )

    # Parser ingesting cluster files
    parser_cluster = subparsers.add_parser(
        "ingest_cluster", help="Indicates that cluster file is being ingested"
    )
    parser_cluster.add_argument(
        "--cluster-file",
        required=True,
        help="Absolute or relative path to cluster file.",
    )
    parser_cluster.add_argument(
        "--ingest-cluster",
        required=True,
        action="store_true",
        help="Indicates that ingest of cluster file should be invoked",
    )
    parser_cluster.add_argument(
        "--name", required=True, help="Name of cluster from input form"
    )
    parser_cluster.add_argument(
        "--domain-ranges",
        type=ast.literal_eval,
        help="Optional paramater taken from UI",
    )

    # Parser ingesting cluster files
    parser_subsample = subparsers.add_parser(
        "ingest_subsample", help="Indicates that subsampling will be initialized"
    )
    parser_subsample.add_argument(
        "--subsample",
        required=True,
        action="store_true",
        help="Indicates that subsampliing functionality should be invoked",
    )
    parser_subsample.add_argument(
        "--cluster-file",
        required=True,
        help="Absolute or relative path to cluster file.",
    )
    parser_subsample.add_argument(
        "--cell-metadata-file", help="Absolute or relative path to cell metadata file."
    )

    return parser


def validate_arguments(parsed_args):
    """Verify parsed input arguments

    Args:
        parsed_args: Parsed input arguments

    Returns:
        None
    """

    if ("matrix_file" in parsed_args and parsed_args.matrix_file_type == "mtx") and (
        parsed_args.gene_file is None or parsed_args.barcode_file is None
    ):
        raise ValueError(
            " Missing arguments: --gene-file and --barcode-file. Mtx files "
            "must include .genes.tsv, and .barcodes.tsv files. See --help for "
            "more information"
        )


def main() -> None:
    """This function handles the actual logic of this script.

    Args:
        None

    Returns:
        None
    """
    status = []
    parsed_args = create_parser().parse_args()
    validate_arguments(parsed_args)
    arguments = vars(parsed_args)
    ingest = IngestPipeline(**arguments)

    if "matrix_file" in arguments:
        ingest.ingest_expression()
    elif "ingest_cell_metadata" in arguments:
        if arguments["ingest_cell_metadata"]:
            status_cell_metadata = ingest.ingest_cell_metadata()
            status.append(status_cell_metadata)
    elif "ingest_cluster" in arguments:
        if arguments["ingest_cluster"]:
            ingest.ingest_cluster()
    elif "subsample" in arguments:
        if arguments["subsample"]:
            ingest.subsample()

    if all(i < 1 for i in status):
        print('ok')
        sys.exit(os.EX_OK)
    else:
        print('bad')
        sys.exit(os.EX_DATAERR)


if __name__ == "__main__":
    main()
