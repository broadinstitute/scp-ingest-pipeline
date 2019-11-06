"""Upload NDJSON file to BigQuery

DESCRIPTION
This CLI takes a local NDJSON file and appends it to an existing bigquery table.
(reference: https://cloud.google.com/bigquery/docs/loading-data-local#loading_data_from_a_local_data_source and
https://cloud.google.com/bigquery/docs/loading-data-cloud-storage-json#loading_json_data_into_a_new_table)

EXAMPLE
$ python upload_to_bq.py dataset_id table_id ../../tests/data/valid_arrays_v1.1.3_for_bq_v1.json

"""

import argparse
from google.cloud import bigquery

client = bigquery.Client()


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention tsv file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('dataset_id', help='<project_ID>:<dataset_id>')
    parser.add_argument('table_id', help='bigquery table_id')
    parser.add_argument('input_json', help='NDJSON file for upload')
    return parser


if __name__ == '__main__':
    args = create_parser().parse_args()
    dataset_id = args.dataset_id
    table_id = args.table_id
    input_json = args.input_json
    dataset_ref = client.dataset(dataset_id)
    table_ref = dataset_ref.table(table_id)
    job_config = bigquery.LoadJobConfig()
    # note WRITE_APPEND is default behavior
    job_config.write_disposition = bigquery.WriteDisposition.WRITE_APPEND
    job_config.source_format = bigquery.SourceFormat.NEWLINE_DELIMITED_JSON
    uri = input_json
    # API request
    with open(input_json, 'rb') as source_file:
        job = client.load_table_from_file(source_file, table_ref, job_config=job_config)

    job.result()  # Waits for table load to complete.
    print("Job finished.")

    print("Loaded {} rows into {}:{}.".format(job.output_rows, dataset_id, table_id))

    destination_table = client.get_table(table_ref)
    print("{} rows in destination table.".format(destination_table.num_rows))
