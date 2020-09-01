"""Helper functions for ingest_pipeline.py
"""

import argparse
import ast

from google.cloud import bigquery
from google.cloud.exceptions import NotFound


# Ingest file types
EXPRESSION_FILE_TYPES = ["dense", "mtx", "loom"]


def bq_dataset_exists(dataset):
    bigquery_client = bigquery.Client()
    dataset_ref = bigquery_client.dataset(dataset)
    exists = False
    try:
        bigquery_client.get_dataset(dataset_ref)
        exists = True
    except NotFound:
        print(f'Dataset {dataset} not found')
    return exists


def bq_table_exists(dataset, table):
    bigquery_client = bigquery.Client()
    dataset_ref = bigquery_client.dataset(dataset)
    table_ref = dataset_ref.table(table)
    exists = False
    try:
        bigquery_client.get_table(table_ref)
        exists = True
    except NotFound:
        print(f'Dataset {table} not found')
    return exists


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
    if "ingest_cell_metadata" in parsed_args:
        if (parsed_args.bq_dataset is not None and parsed_args.bq_table is None) or (
            parsed_args.bq_dataset is None and parsed_args.bq_table is not None
        ):
            raise ValueError(
                'Missing argument: --bq_dataset and --bq_table are both required for BigQuery upload.'
            )
        if parsed_args.bq_dataset is not None and not bq_dataset_exists(
            parsed_args.bq_dataset
        ):
            raise ValueError(
                f' Invalid argument: unable to connect to a BigQuery dataset called {parsed_args.bq_dataset}.'
            )
        if parsed_args.bq_table is not None and not bq_table_exists(
            parsed_args.bq_dataset, parsed_args.bq_table
        ):
            raise ValueError(
                f' Invalid argument: unable to connect to a BigQuery table called {parsed_args.bq_table}.'
            )


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
        "--study-file-id",
        required=True,
        help="Single study accession associated with ingest files",
    )
    parser.add_argument("--study-id", required=True, help="MongoDB identifier")

    profile_memory_text = "Whether to profile memory.  Outputs \
        mprofile_<YYYYMMDDHHMMSS>.dat file locally.  \
        Use locally for optimization work. \
        See https://github.com/pythonprofilers/memory_profiler#time-based-memory-usage"
    parser.add_argument(
        "--profile-memory", action="store_true", help=profile_memory_text
    )

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
        "--study-accession",
        required=True,
        help="Single study accession associated with ingest files.",
    )
    parser_cell_metadata.add_argument(
        "--bq-dataset", help="BigQuery dataset identifer for ingest job."
    )
    parser_cell_metadata.add_argument(
        "--bq-table", help="BigQuery table identifer for ingest job."
    )
    parser_cell_metadata.add_argument(
        "--ingest-cell-metadata",
        required=True,
        action="store_true",
        help="Indicates that ingest of cell metadata should be invoked",
    )
    parser_cell_metadata.add_argument(
        "--validate-convention",
        action="store_true",
        help="Indicates that metadata file should be validated against convention",
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
        "--name", required=True, help="Name of cluster from input form"
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
