"""Launch ingest pipeline stress test

DESCRIPTION
(Internal use only) Launch single study creation followed by file upload using manage-study
OR concurrent ingest requests for N studies (stress test)

PREREQUISITES
set up access token on command line (Script requires access to non-public Google bucket)
ACCESS_TOKEN=`gcloud auth print-access-token`

EXAMPLES
Creates N new study/ies in staging
python ingest_stress_test.py --token=$ACCESS_TOKEN -n <number of studies>

Uses existing study, 'stress-devel-study' (SCP138 in staging), uploads 8K data set
# Usage conditions: must have VPN access to staging server
#                   files being uploaded cannot already exist in createTest
python ingest_stress_test.py --token=$ACCESS_TOKEN --dev-run
"""

import argparse
import subprocess
from datetime import datetime

DATA = {
    # matrix has ~8K cells, ~20K genes
    # baseline ingest for dense 3.6m
    "small": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/small_expr.txt.gz",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/small_coords.txt",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/small_metadata.txt",
    },
    # matrix has ~18.7K cells, ~25K genes
    # baseline ingest for dense 4.5m
    "medium": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/medium_expr.txt.gz",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/medium_coords.txt",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/medium_metadata.txt",
    },
    # matrix has ~220K cells, ~16K genes
    # baseline ingest for dense ~25.5m
    "large": {
        "dense": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/large_expr.txt",
        "cluster": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/large_coords.txt",
        "metadata": "gs://fc-2f8ef4c0-b7eb-44b1-96fe-a07f0ea9a982/test_Data/scp-ingest-pipeline/stress_test_data/large_metadata.txt",
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
        help='Number of studies to create for stress test',
        type=int,
        default=1,
    )

    parser.add_argument(
        '--environment',
        choices=['development', 'staging', 'production'],
        help='Server to test against',
        default='staging',
    )

    parser.add_argument(
        '--token',
        default=None,
        required=True,
        help='Personal token after logging into Google (OAuth2).  This token is not persisted after the finish of the script.',
    )

    parser.add_argument(
        '--dev-run',
        action='store_true',
        help='Run using existing study, "stress-devel-study" (SCP138)',
    )

    parser.add_argument(
        '--data-set',
        choices=['small', 'medium', 'large'],
        default='small',
        help='Data set for stress test',
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help='Show upload command results for troubleshooting',
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show upload commands without running upload request',
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
                'test study for smoke or stress test',
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
                "GRCh37",
            ],
        }

        if args.dry_run:

            class Result:
                def __init__(self):
                    self.returncode = 0
                    self.stdout = b''

            result = Result()
            print(f'{self.action_setup} {study_actions.get(action)}')
        else:
            result = subprocess.run(
                self.action_setup + study_actions.get(action),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        if args.debug:
            print(result)
        if result.returncode == 0:
            print(f"Success {result.stdout.decode('utf-8')}")
            return True
        # for "study_exists" action, failure to confirm a study exists is not an error
        elif action == "study_exists":
            return False
        elif result.stdout.decode('utf-8') == "Error 403: Forbidden\n":
            print(result.stdout.decode('utf-8'))
            print("Check that you have access (ie. VPN) to the Portal server")
            return False
        else:
            print(
                f"Error :{result.stdout.decode('utf-8')} {result.stderr.decode('utf-8')}"
            )
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
    study_names = [f'{test_type}-{test_id}-{i}' for i in range(args.number)]

    # create studies unless it is a dev-run
    if not args.dev_run and not msaction.send_request(study_names[0], "study_exists"):
        for study in study_names:
            print(f'Requesting creation of {study}')
            create = msaction.send_request(study, "create_study")
            if not create:
                print(f'Failed to create {study}. Exiting')
                exit(1)

    # if a dev-run, limit ingest jobs to single, pre-created test study
    if args.dev_run:
        study_names = ['stress-devel-study']

    print(f'Using the following study/ies for file upload {study_names}')
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
            print(f'Expression matrix upload failed for {study}. Exiting')
            exit(1)
        upload_cluster = msaction.send_request(study, "upload_cluster")
        if not upload_cluster:
            print(f'Cluster upload failed for {study}. Exiting')
            exit(1)
        upload_metadata = msaction.send_request(study, "upload_metadata")
        if not upload_metadata:
            print(f'Metadata upload failed for {study}. Exiting')
            exit(1)
