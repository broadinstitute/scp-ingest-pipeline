# scp-ingest-pipeline

File Ingest Pipeline for Single Cell Portal

[![Build status](https://img.shields.io/circleci/build/github/broadinstitute/scp-ingest-pipeline.svg)](https://circleci.com/gh/broadinstitute/scp-ingest-pipeline)
[![Code coverage](https://codecov.io/gh/broadinstitute/scp-ingest-pipeline/branch/master/graph/badge.svg)](https://codecov.io/gh/broadinstitute/scp-ingest-pipeline)

The SCP Ingest Pipeline is an ETL pipeline for single-cell RNA-seq data.

# Prerequisites

- Python 3.10
- Google Cloud Platform project
- Suitable service account (SA) and MongoDB VM in GCP. SA needs roles "Editor", "Genomics Pipelines Runner", and "Storage Object Admin". Broad Institute engineers: see instructions [here](https://github.com/broadinstitute/single_cell_portal_configs/tree/master/terraform-mongodb).
- SAMtools, if using `ingest/make_toy_data.py`
- Tabix, if using `ingest/genomes/genomes_pipeline.py`

# Install

### Native

Fetch the code, boot your virtualenv, install dependencies:

```
git clone git@github.com:broadinstitute/scp-ingest-pipeline.git
cd scp-ingest-pipeline
python3 -m venv env --copies
source env/bin/activate
pip install -r requirements.txt
source scripts/setup-mongo-dev.sh
```

### Docker

With Docker running and `gcloud` authenticated on your local machine, run:

```
scripts/docker-compose-setup.sh
```

If on Apple silicon Mac (e.g. M1), and performance seems poor, consider generating a docker image using the arm64 base. Example test image: gcr.io/broad-singlecellportal-staging/single-cell-portal:development-2.2.0-arm64, usage:

```
scripts/docker-compose-setup.sh -i development-2.2.0-arm64
```

To update dependencies when in Docker, you can pip install from within the Docker Bash shell after adjusting your requirements.txt.
If you close your shell after that, your newly installed dependencies will be lost.  Dependencies only persist after merging your
new requirements.txt into `development`.  TODO (SCP-4941): Add entry-point script to run `pip install`.

### Optional

To use `ingest/make_toy_data.py`:

```
brew install samtools
```

To use `ingest/genomes/genomes_pipeline.py`:

```
brew install tabix
```

## Git hooks

After installing Ingest Pipeline, add Git hooks to help ensure code quality:

```
pre-commit install && pre-commit install -t pre-push
```

The hooks will expect that [git-secrets](https://github.com/awslabs/git-secrets) has been set up. If you are a Broad Institute employee who has not done this yet, please see: [broadinstitute/single_cell_portal_configs](https://github.com/broadinstitute/single_cell_portal_configs) for specific guidance.

### Bypass hooks

In rare cases, you might need to skip Git hooks, like so:

- Skip commit hooks: `git commit ... --no-verify`
- Skip pre-push hooks: `git push ... --no-verify`

# Test

After [installing](#install):

```
source env/bin/activate
cd tests

# Run all tests
pytest
```

Some common `pytest` usage examples (run in `/tests`):

```
# Run all tests and see print() output
pytest -s

# Run only tests in test_ingest.py
pytest test_ingest.py

# Run all tests, show code coverage metrics
pytest --cov=../ingest/
```

For more, see <https://docs.pytest.org/en/stable/usage.html>.

## Testing in Docker locally
<!--
Step 1 is also useful for troubleshooting when Dockerfile updates fail to build
-->
If you have difficulties installing and configuring `scp-ingest-pipeline` due to hardware issues (e.g. Mac M1 chips),
you can alternatively test locally by building the Docker image and then running any commands inside the container.
There are some extra steps required, but this sidesteps the need to install packages locally.

### 1. Build the image

Run the following command to build the testing Docker image locally (make sure Docker is running first). This build command will incorporate any changes in the local instance of your repo, committed or not:

```
docker build -t ingest-pipeline:test-candidate .
```

Note - the base Ubuntu image used in the Dockerfile comes from Google Container Registry (GCR, aka gcr.io), if this is your first time doing `docker build` you may need to configure Docker to use the Google Cloud CLI to authenticate requests to Container Registry:

```
gcloud auth configure-docker
```

Pro-Tip: For local builds, you can try adding docker build options `--progress=plain` (for more verbose build info) and/or `--no-cache` (when you want to ensure a build with NO cached layers)

### 2. Set up environment variables

Run the following to pull database-specific secrets out of Google Secrets Manager (GSM):

```
source scripts/setup-mongo-dev.sh
```

Now run `env` to make sure you've set the following values:

```
MONGODB_USERNAME=single_cell
DATABASE_NAME=single_cell_portal_development
MONGODB_PASSWORD=<password>
DATABASE_HOST=<ip address>
```

### 3. Print out your service account keyfile

Run the following to export out your default service account JSON keyfile:

```
GOOGLE_PROJECT=$(gcloud info --format="value(config.project)")
gcloud secrets versions access latest --project=$GOOGLE_PROJECT --secret=default-sa-keyfile | jq > /tmp/keyfile.json
```

### 4. Start the Docker container

Run the container, passing in the proper environment variables:

```
docker run --name scp-ingest-test -e MONGODB_USERNAME="$MONGODB_USERNAME" -e DATABASE_NAME="$DATABASE_NAME" \
           -e MONGODB_PASSWORD="$MONGODB_PASSWORD" -e DATABASE_HOST="$DATABASE_HOST" \
           -e GOOGLE_APPLICATION_CREDENTIALS=/tmp/keyfile.json --rm -it \
           ingest-pipeline:test-candidate bash
```

Note: on an M1 machine, you may see this message:

> WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested

### 5. Copy keyfile to running container

In a separate terminal window, copy the JSON keyfile from above to the expected location:

```
docker cp /tmp/keyfile.json scp-ingest-test:/tmp
```

You can now run any `ingest_pipeline.py` command you wish inside the container.

# Use

Run this every time you start a new terminal to work on this project:

```
source env/bin/activate
```

See [`ingest_pipeline.py`](https://github.com/broadinstitute/scp-ingest-pipeline/blob/development/ingest/ingest_pipeline.py) for usage examples.

## Troubleshooting during set up

If you run into an error like: "... [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed ... " try:

- Open terminal
- `cd` to where python is installed
- Run the certificates command with `/Applications/Python\ < Your Version of Python Here >/Install\ Certificates.command`

If you run into an error like "ModuleNotFoundError: No module named 'google'" try:

- Open terminal
- Run `pip install --upgrade google-api-python-client`
