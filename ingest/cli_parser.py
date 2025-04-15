"""Helper functions for ingest_pipeline.py"""

import argparse
import ast

from google.cloud import bigquery
from google.cloud.exceptions import NotFound


# Ingest file types
EXPRESSION_FILE_TYPES = ["dense", "mtx", "h5ad"]


def bq_dataset_exists(dataset):
    bigquery_client = bigquery.Client()
    dataset_ref = bigquery_client.dataset(dataset)
    exists = False
    try:
        bigquery_client.get_dataset(dataset_ref)
        exists = True
    except NotFound:
        print(f"Dataset {dataset} not found")
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
        print(f"Dataset {table} not found")
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
                "Missing argument: --bq_dataset and --bq_table are both required for BigQuery upload."
            )
        if parsed_args.bq_dataset is not None and not bq_dataset_exists(
            parsed_args.bq_dataset
        ):
            raise ValueError(
                f" Invalid argument: unable to connect to a BigQuery dataset called {parsed_args.bq_dataset}."
            )
        if parsed_args.bq_table is not None and not bq_table_exists(
            parsed_args.bq_dataset, parsed_args.bq_table
        ):
            raise ValueError(
                f" Invalid argument: unable to connect to a BigQuery table called {parsed_args.bq_table}."
            )
    if (
        "differential_expression" in parsed_args
        and parsed_args.annotation_type != "group"
    ):
        raise ValueError(
            "Differential expression analysis restricted to group-type annotations,"
            f" cannot run on data of type \"{parsed_args.annotation_type}\"."
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

    parser.add_argument(
        "--user-metrics-uuid",
        required=False,
        type=is_valid_uuid,
        help="User identifier for Bard, i.e. the user's Mixpanel distinct ID",
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
    parser_cell_metadata.add_argument(
        "--has-modality",
        type=ast.literal_eval,
        help="Array of modalities metadata ingest should transform to boolean for BigQuery",
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

    # Differential expression subparsers
    parser_differential_expression = subparsers.add_parser(
        "differential_expression",
        help="Indicates differential expression analysis processing",
    )

    parser_differential_expression.add_argument(
        "--differential-expression",
        required=True,
        action="store_true",
        help="Indicates that differential expression analysis should be invoked",
    )

    parser_differential_expression.add_argument(
        "--de-type",
        default="rest",
        choices=['rest', 'pairwise'],
        help="Accepted values: 'pairwise' or 'rest' (default)",
    )

    parser_differential_expression.add_argument(
        "--raw-location",
        help="location of raw counts. '.raw' for raw slot, "
        "else adata.layers key value",
    )

    parser_differential_expression.add_argument(
        "--study-accession",
        required=True,
        help="Single study accession associated with provided DE input files.",
    )

    parser_differential_expression.add_argument(
        "--annotation-name", required=True, help="Name of annotation for DE analysis"
    )

    parser_differential_expression.add_argument(
        "--annotation-type", required=True, help="Type of annotation for DE analysis"
    )

    parser_differential_expression.add_argument(
        "--annotation-scope", required=True, help="Scope of annotation for DE analysis"
    )

    parser_differential_expression.add_argument(
        "--method", default="wilcoxon", help="method for DE"
    )

    parser_differential_expression.add_argument(
        "--cluster-name", required=True, help="study owner-specified cluster name"
    )

    parser_differential_expression.add_argument(
        "--cluster-file",
        required=True,
        help="Absolute or relative path to cluster file.",
    )

    parser_differential_expression.add_argument(
        "--annotation-file",
        required=True,
        help="Absolute or relative path to cell metadata or cluster file of annotations.",
    )

    parser_differential_expression.add_argument(
        "--matrix-file-path",
        required=True,
        help="Absolute or relative path to "
        "expression file. For 10x data this is "
        "the .mtx file",
    )

    parser_differential_expression.add_argument(
        "--matrix-file-type",
        choices=EXPRESSION_FILE_TYPES,
        type=str.lower,
        required=True,
        help=matrix_file_type_txt,
    )

    # Gene and Barcode arguments for MTX bundle
    parser_differential_expression.add_argument(
        "--barcode-file", help="Path to .barcodes.tsv files"
    )
    parser_differential_expression.add_argument(
        "--gene-file", help="Path to .genes.tsv file"
    )

    # For pairwise analyses
    parser_differential_expression.add_argument(
        "--group1",
        help="1st annotation label to use for pairwise DE analysis",
    )

    parser_differential_expression.add_argument(
        "--group2",
        help="2nd annotation label to use for pairwise DE analysis",
    )

    parser_ingest_differential_expression = subparsers.add_parser(
        "ingest_differential_expression",
        help="Indicates author differential expression analysis processing",
    )

    parser_ingest_differential_expression.add_argument(
        "--ingest-differential-expression",
        required=True,
        action="store_true",
        help="Indicates that author differential expression analysis should be invoked",
    )

    parser_ingest_differential_expression.add_argument(
        "--study-accession",
        required=True,
        help="Single study accession associated with provided DE input files.",
    )

    parser_ingest_differential_expression.add_argument(
        "--annotation-name", required=True, help="Name of annotation for DE analysis"
    )

    parser_ingest_differential_expression.add_argument(
        "--annotation-type", required=True, help="Type of annotation for DE analysis"
    )

    parser_ingest_differential_expression.add_argument(
        "--annotation-scope", required=True, help="Scope of annotation for DE analysis"
    )

    parser_ingest_differential_expression.add_argument(
        "--method", default="wilcoxon", help="method for DE"
    )

    parser_ingest_differential_expression.add_argument(
        "--cluster-name", required=True, help="study owner-specified cluster name"
    )

    parser_ingest_differential_expression.add_argument(
        "--differential-expression-file",
        required=True,
        help="Path to DE file uploaded by author.",
    )

    parser_ingest_differential_expression.add_argument(
        "--gene-header",
        required=True,
        help="Header used for gene names / symbols in DE file",
    )

    parser_ingest_differential_expression.add_argument(
        "--group-header", required=True, help="Header used for group in DE file"
    )

    parser_ingest_differential_expression.add_argument(
        "--comparison-group-header",
        required=False,
        help=(
            "Header used for comparison group in DE file.  "
            + "For pairwise comparisons.  Can omit if DE file is in one-vs-rest-only format."
        ),
    )

    parser_ingest_differential_expression.add_argument(
        "--size-metric",
        required=True,
        help='Header used as size metric in DE file, e.g. "logfoldchanges", "avg_log2FC", etc.',
    )

    parser_ingest_differential_expression.add_argument(
        "--significance-metric",
        required=True,
        help='Header used as significance metric in DE file, e.g. "pvals_adj", "p_val_adj", etc.',
    )

    # AnnData subparsers
    parser_anndata = subparsers.add_parser(
        "ingest_anndata", help="Indicates that AnnData file is being ingested"
    )

    parser_anndata.add_argument(
        "--ingest-anndata",
        required=True,
        action="store_true",
        help="Indicates that ingest of AnnData file should be invoked",
    )

    parser_anndata.add_argument(
        "--anndata-file", required=True, help="Path to AnnData file"
    )

    parser_anndata.add_argument(
        "--obsm-keys",
        type=ast.literal_eval,
        help="Array of obsm key(s) to extract as cluster files",
    )

    parser_anndata.add_argument(
        "--raw-location",
        help="location of raw counts. '.raw' for raw slot, "
        "else adata.layers key value or None if no raw counts",
    )

    parser_anndata.add_argument(
        "--extract",
        type=ast.literal_eval,
        help="Array of file types to extract, options include ['cluster', 'metadata', 'processed_expression']",
    )

    parser_expression_writer = subparsers.add_parser(
        "render_expression_arrays",
        help="Indicates preprocessing of cluster/expression files for image pipeline",
    )

    parser_expression_writer.add_argument(
        '--render-expression-arrays',
        action="store_true",
        help='Invoke expression_writer.py',
        required=True,
    )

    parser_expression_writer.add_argument(
        '--cluster-file', help='path to cluster file', required=True
    )
    parser_expression_writer.add_argument(
        '--cluster-name', help='name of cluster object', required=True
    )
    parser_expression_writer.add_argument(
        '--matrix-file-path', help='path to matrix file', required=True
    )
    parser_expression_writer.add_argument(
        '--matrix-file-type',
        help='type to matrix file (dense or mtx)',
        required=True,
        choices=['dense', 'mtx'],
    )
    parser_expression_writer.add_argument(
        '--gene-file', help='path to gene file (omit for dense matrix files)'
    )
    parser_expression_writer.add_argument(
        '--barcode-file', help='path to barcode file (omit for dense matrix files)'
    )

    parser_rank_genes = subparsers.add_parser(
        "rank_genes",
        help="Rank genes in a study by mentions in publication, DE, and global interest",
    )

    parser_rank_genes.add_argument(
        '--rank-genes',
        action="store_true",
        help='Invoke rank_genes.py',
        required=True,
    )

    parser_rank_genes.add_argument(
        '--study-accession', help='Study accession, e.g. "SCP123"', required=True
    )
    parser_rank_genes.add_argument(
        '--publication',
        help="URL of the study's publicly-accessible research article, or GS URL or local path to publication text file",
        required=True,
    )

    return parser


def is_valid_uuid(value):
    import uuid

    try:
        if value:
            uuid.UUID(value)
            return value
        else:
            return None
    except ValueError as e:
        raise argparse.ArgumentTypeError(e)
