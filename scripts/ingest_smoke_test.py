"""Launch ingest pipeline smoke test

DESCRIPTION
Launch single study creation followed by file upload using manage-study
OR concurrent ingest requests for X studies (stress test)

PREREQUISITES
set up access token on command line
ACCESS_TOKEN=`gcloud auth print-access-token`

EXAMPLE
Creates N new study/ies in staging
python ingest_smoke_test.py --token=$ACCESS_TOKEN -n <number of studies>

Uses existing study, 'createTest' (SCP88 in staging), uploads 8K data set
# Usage conditions: must have VPN access to staging server
#                   files being uploaded cannot already exist in createTest
python ingest_smoke_test.py --token=$ACCESS_TOKEN --dev-run
"""

import argparse
import subprocess
from datetime import datetime

DATA = {
    # matrix has ~8K cells, ~20K genes (aka "mult-bcl" data set); SCP23 and SCP31 in staging
    # baseline ingest for dense TBD
    "small": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/mult_bcl_out/mult_bcl_out.scp.expr.txt.gz",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/mult_bcl_out/mult_bcl_out.scp.X_tsne.coords.txt",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/mult_bcl_out/mult_bcl_out.scp.metadata.txt",
    },
    # matrix has ~18.7K cells, ~25K genes (aka "PolypScrape" data from Nasal polyps)
    # baseline ingest for dense ~4.5m
    "medium": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/test_data_for_production_ingest/nasal_polyps/20180921_PolypScrape_cleaned_data.txt.gz",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/test_data_for_production_ingest/nasal_polyps/20180921_PolypScrape_cleaned_tsnecoord.txt",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/test_data_for_production_ingest/nasal_polyps/alexandria_structured_metadata.txt",
    },
    # matrix has ~220K cells, ~16K genes from SCP947 in prod
    # baseline ingest for dense ~25.5m
    "large": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/SCP947/All.exp.tsv",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/SCP947/All.cluster.tsv",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/launch_ingest_test_data/SCP947/metadata.tsv",
    },
}


def create_parser():
    """Parse command line values
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-n',
        '--number',
        help='number of studies to create for stress test',
        type=int,
        default=1,
    )

    parser.add_argument(
        '--environment',
        help='server to test against [development,staging,production]',
        default='staging',
    )

    parser.add_argument(
        '--token',
        default=None,
        help='Personal token after logging into Google (OAuth2).  This token is not persisted after the finish of the script.',
    )

    parser.add_argument(
        '--dev-run',
        action='store_true',
        help='Run using stress_20200810-1357030 as study.',
    )

    parser.add_argument(
        '--data-set',
        default='small',
        help='Select data set for stress test [small, medium, large].',
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show upload commands without running upload request.',
    )

    return parser


class ManageStudyAction:
    def __init__(self, data_set):
        self.data_set = data_set
        self.action_setup = [
            'manage-study',
            '--environment',
            args.environment,
            '--token',
            args.token,
        ]

    def send_request(self, study, action):

        study_actions = {
            "study_exists": [
                'get-study-attribute',
                '--study-name',
                study,
                '--attribute',
                'accession',
            ],
            'create_study': [
                'create-study',
                '--study-name',
                study,
                '--description',
                'stress test study',
                '--is-private',
            ],
            "upload_cluster": [
                'upload-cluster',
                '--study-name',
                study,
                '--file',
                DATA[self.data_set]['cluster'],
                '--cluster-name',
                'tSNE',
            ],
            "upload_metadata": [
                'upload-metadata',
                '--study-name',
                study,
                '--file',
                DATA[self.data_set]['metadata'],
            ],
            "upload_dense_expression": [
                'upload-expression',
                '--study-name',
                study,
                '--file',
                DATA[self.data_set]['dense'],
                '--species',
                "human",
                '--genome',
                "hg19",
            ],
        }

        if args.dry_run:

            class Result:
                def __init__(self):
                    self.returncode = 0
                    self.stdout = b''

            result = Result()
            print(self.action_setup + study_actions.get(action))
        else:
            result = subprocess.run(
                self.action_setup + study_actions.get(action),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

        print(result)
        if result.returncode == 0:
            print('Success', result.stdout.decode('utf-8'))
            return True
        # for "study_exists", failure to confirm a study exists is not an error
        elif action == "study_exists":
            return False
        else:
            print('Error', result.stdout.decode('utf-8'), result.stderr.decode('utf-8'))
            return False


if __name__ == '__main__':
    args = create_parser().parse_args()

    msaction = ManageStudyAction(args.data_set)

    # generate study names
    test_id = datetime.now().strftime('%Y%m%d-%H%M%S')
    if args.number == 1:
        test_type = "smoke"
    else:
        test_type = "stress"
    study_names = [test_type + '-' + test_id + '-' + str(i) for i in range(args.number)]

    # create studies unless it is a dev-run
    if not args.dev_run and not msaction.send_request(study_names[0], "study_exists"):
        for study in study_names:
            print('Requesting creation of', study)
            create = msaction.send_request(study, "create_study")
            if not create:
                print('Failed to create', study.Exiting)
                exit(1)

    # if a dev-run, limit ingest jobs to single, pre-created test study
    if args.dev_run:
        study_names = ['createTest']

    print('Using the following study/ies for file upload', study_names)
    # obtain study file info

    # launch ingest jobs for expression matrix, metadata and cluster files
    # exits upon first ingest request failure
    print('Requesting file upload for')
    for study in study_names:
        print(study)
        upload_dense_expression = msaction.send_request(
            study, "upload_dense_expression"
        )
        if not upload_dense_expression:
            print('Expression upload failed for ', study, ". Exiting")
            exit(1)
        upload_cluster = msaction.send_request(study, "upload_cluster")
        if not upload_cluster:
            print('Cluster upload failed for ', study, ". Exiting")
            exit(1)
        upload_metadata = msaction.send_request(study, "upload_metadata")
        if not upload_metadata:
            print('Metadata upload failed for ', study, ". Exiting")
            exit(1)

