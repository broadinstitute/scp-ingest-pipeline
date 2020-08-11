"""Pipeline to download, transform, and upload genome reference data for SCP

SUMMARY

This ETL pipeline downloads, transforms, and uploads reference data for
genome assemblies and genome annotations needed for Single Cell Portal (SCP).

Raw data is fetched from NCBI and Ensembl, processed to transform it as needed,
then uploaded to a Google Cloud Storage (GCS) bucket for SCP.  The uploaded data
is used to display interactive genome visualizations, e.g. igv.js and
Ideogram.js, and (in the future) to run biological workflows / analysis
pipelines, e.g. Cell Ranger and inferCNV.

To use this reference data in SCP, upload the "species_reference_metadata.tsv"
file in the "outputs" directory using the "Upload Species List" button in
https://portals.broadinstitute.org/single_cell/species.

INSTALL

python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt

cd ingest/genomes

EXAMPLES

# Basic usage.  Upload GTF products to reference_dev in default SCP GCS bucket.
$ python3 genomes_pipeline.py --vault-path=secrets/service_account.json

# Upload GTF products to reference_data_staging folder, using cached data from previous run
$ python3 genomes_pipeline.py --vault-path=secrets/service_account.json --input-dir ../../../lib/assets/python/genomes/ --remote-output-dir reference_data_staging/ --use-cache
"""

import argparse

from genome_annotations import GenomeAnnotations
from genome_assemblies import GenomeAssemblies


def create_parser() -> argparse.ArgumentParser:
    """Creates parser for input arguments.

    Structuring the argument parsing code like this eases automated testing.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--vault-path', help='Path in Vault for GCS service account credentials'
    )
    parser.add_argument(
        '--input-dir',
        help='Input directory; where to find organisms.tsv.  Default: ./',
        default='./',
    )
    parser.add_argument(
        '--remote-output-dir',
        help='Remote directory for output in GCS bucket.  '
        + 'Default: reference_data_dev/',
        default='reference_data_dev/',
    )
    parser.add_argument(
        '--local-output-dir',
        help='Local directory for output.  Default: output/',
        default='output/',
    )
    parser.add_argument(
        '--gcs-bucket',
        help='Name of GCS bucket for upload.  ' + 'Default: reference_data_dev/',
        default='fc-bcc55e6c-bec3-4b2e-9fb2-5e1526ddfcd2',
    )
    parser.add_argument(
        '--copy-data-from-prod-dir',
        help='Remote directory from which to copy data into '
        + 'remote_output_dir.  Use to ensure test data '
        + 'environment is equivalent to production data '
        + 'environment.  Default: reference_data/',
        default='reference_data/',
    )
    parser.add_argument('--use-cache', help='Whether to use cache', action='store_true')

    return parser


def parse_assemblies(args):
    GenomeAssemblies(input_dir=args.input_dir, output_dir=args.local_output_dir)


def parse_genome_annotations(args):
    """Download genome annotations, transform for visualizations, upload to GCS
    """
    input_dir = args.input_dir
    remote_output_dir = args.remote_output_dir
    vault_path = args.vault_path
    local_output_dir = args.local_output_dir
    gcs_bucket = args.gcs_bucket
    remote_output_dir = args.copy_data_from_prod_dir
    use_cache = args.use_cache

    GenomeAnnotations(
        vault_path,
        input_dir,
        local_output_dir,
        gcs_bucket,
        remote_output_dir,
        use_cache,
    ).load()


def main() -> None:
    """Enables running via module or CLI
    """
    args = create_parser().parse_args()
    parse_assemblies(args)
    parse_genome_annotations(args)


if __name__ == "__main__":
    main()
