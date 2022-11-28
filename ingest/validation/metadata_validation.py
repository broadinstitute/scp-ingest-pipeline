"""Validate input metadata TSV file against metadata convention.

DESCRIPTION
This script takes a TSV metadata file and validates against a metadata convention
using the python jsonschema library. The metadata convention JSON schema
represents the rules that should be enforced on metadata files for studies
participating under the convention.

EXAMPLE
# Using JSON file for latest Alexandria metadata convention in repo, validate input TSV
$ python3 metadata_validation.py  ../../tests/data/valid_no_array_v2.0.0.tsv

# generate an issues.json file to compare with reference test files
$ python3 metadata_validation.py --issues-json ../../tests/data/valid_no_array_v2.0.0.tsv

# generate a BigQuery upload file to compare with reference test files
$ python3 metadata_validation.py --bq-json ../../tests/data/valid_no_array_v2.0.0.tsv

# use a different metadata convention for validation
$ python3 metadata_validation.py --convention <path to convention json> ../../tests/data/valid_no_array_v2.0.0.tsv

"""

import argparse
import json
import logging
from collections import defaultdict
import sys
import requests
import urllib.parse as encoder
import re
import os
import numbers
import time
import backoff
import csv
import copy
import itertools
import math
import pandas as pd

import colorama
from colorama import Fore
import jsonschema
from google.cloud import bigquery

sys.path.append("..")
try:
    # Used when importing internally and in tests
    from cell_metadata import CellMetadata
    from validation.validate_metadata import (
        report_issues,
        validate_input_metadata,
        write_metadata_to_bq,
        serialize_issues,
        exit_if_errors,
    )
except ImportError:
    # Used when importing as external package, e.g. imports in single_cell_portal code
    from ..cell_metadata import CellMetadata


def create_parser():
    """Parse command line values for validate_metadata
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # parser.add_argument(
    #     '--output',
    #     '-o',
    #     type=str,
    #     help='Output file name [optional]',
    #     default=None
    # )
    # parser.add_argument(
    #     '--key_id',
    #     '-k'
    #     type=str,
    #     help='Key metadata name for parsing; CellID for metadata, BiosampleID for sample sheets [optional]',
    #     default='CellID'
    # )

    # helper param to create JSON representation of metadata.issues
    # to generate reference output for tests
    parser.add_argument("--issues-json", action="store_true")
    # helper param to create JSON representation of convention metadata
    # to generate json for bigquery testing
    parser.add_argument("--bq-json", action="store_true")
    # overwrite existing output
    parser.add_argument("--force", action="store_true")
    # test BigQuery upload functions
    parser.add_argument("--upload", action="store_true")
    # validate_metadata.py CLI only for dev, bogus defaults below shouldn't propagate
    # make bogus defaults obviously artificial for ease of detection
    parser.add_argument(
        "--study-id",
        help="MongoDB study identifier",
        default="dec0dedfeed1111111111111",
    )
    parser.add_argument(
        "--study-file-id",
        help="MongoDB file identifier",
        default="addedfeed000000000000000",
    )
    parser.add_argument(
        "--study-accession", help="SCP study accession", default="SCPtest"
    )
    parser.add_argument(
        "--bq-dataset", help="BigQuery dataset identifier", default="cell_metadata"
    )
    parser.add_argument(
        "--bq-table", help="BigQuery table identifier", default="alexandria_convention"
    )
    parser.add_argument(
        "--convention",
        help="Metadata convention JSON file",
        default="../../schema/alexandria_convention/alexandria_convention_schema.json",
    )
    parser.add_argument("input_metadata", help="Metadata TSV file")
    return parser


def check_if_old_output():
    """Exit if old output files found
    """
    output_files = ["bq.json"]

    old_output = False
    for file in output_files:
        if os.path.exists(file):
            print(f"{file} already exists, please delete file and try again")
            old_output = True
    if old_output:
        exit(1)


if __name__ == "__main__":
    args = create_parser().parse_args()
    arguments = vars(args)
    if not args.force:
        check_if_old_output()

    with open(args.convention, "r") as f:
        convention = json.load(f)
    metadata = CellMetadata(
        file_path=args.input_metadata,
        study_id=args.study_id,
        study_file_id=args.study_file_id,
        study_accession=args.study_accession,
    )
    metadata.preprocess(True)
    print("Validating", args.input_metadata)

    validate_input_metadata(metadata, convention, args.bq_json)
    if args.issues_json:
        serialize_issues(metadata)
    report_issues(metadata)
    if args.upload:
        write_metadata_to_bq(metadata, args.bq_dataset, args.bq_table)
    exit_if_errors(metadata)
