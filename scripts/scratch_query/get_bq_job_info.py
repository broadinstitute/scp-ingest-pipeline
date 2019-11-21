"""Upload NDJSON file to BigQuery

DESCRIPTION
This CLI takes a local NDJSON file and appends it to an existing bigquery table.
(reference: https://cloud.google.com/bigquery/docs/managing-jobs
https://cloud.google.com/bigquery/troubleshooting-errors)

SYNTAX
$ python query_bq_job.py job_id

EXAMPLE
$ python query_bq_job.py d30e45f5-da00-4dc5-8c52-539f67274725

"""

import argparse
from google.cloud import bigquery


def create_parser():
    """
    Command Line parser for serialize_convention

    Input: metadata convention tsv file
    """
    # create the argument parser
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('job_id', help='bigqury job ID')
    return parser


if __name__ == '__main__':
    args = create_parser().parse_args()
    job_id = args.job_id
    client = bigquery.Client()

    location = 'us'

    job = client.get_job(job_id, location=location)

    print("Details for job {} running in {}:".format(job_id, location))
    print(
        "\tType: {}\n\tState: {}\n\tCreated: {}".format(
            job.job_type, job.state, job.created, job.status
        )
    )
