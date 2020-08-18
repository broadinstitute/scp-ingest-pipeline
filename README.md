# scp-ingest-pipeline

File Ingest Pipeline for Single Cell Portal

[![Build status](https://img.shields.io/circleci/build/github/broadinstitute/scp-ingest-pipeline.svg)](https://circleci.com/gh/broadinstitute/scp-ingest-pipeline)
[![Code coverage](https://codecov.io/gh/broadinstitute/scp-ingest-pipeline/branch/master/graph/badge.svg)](https://codecov.io/gh/broadinstitute/scp-ingest-pipeline)

The SCP Ingest Pipeline is an ETL pipeline for single-cell RNA-seq data.

# Prerequisites

- Python 3.7+
- Google Cloud Platform project
- Suitable service account (SA) and MongoDB VM in GCP. SA needs roles "Editor", "Genomics Pipelines Runner", and "Storage Object Admin". Broad Institute engineers: see instructions [here](https://github.com/broadinstitute/single_cell_portal_configs/tree/master/terraform-mongodb).
- SAMtools, if using `ingest/make_toy_data.py`
- Tabix, if using `ingest/genomes/genomes_pipeline.py`

# Install

Fetch the code, boot your virtualenv, install dependencies:

```
git clone git@github.com:broadinstitute/scp-ingest-pipeline.git
cd scp-ingest-pipeline
python3 -m venv env --copies
source env/bin/activate
pip install -r requirements.txt
```

To use `ingest/make_toy_data.py`:

```
brew install samtools
```

To use `ingest/genomes/genomes_pipeline.py`:

```
brew install tabix
```

Now get secrets from Vault to set environment variables needed to write to the database:

```
export BROAD_USER="<username in your email address>"

export DATABASE_NAME="single_cell_portal_development"

vault login -method=github token=`~/bin/git-vault-token`

# Get username and password
vault read secret/kdux/scp/development/$BROAD_USER/mongo/user

export MONGODB_USERNAME="<username from Vault>"
export MONGODB_PASSWORD="<password from Vault>"

# Get external IP address for host
vault read secret/kdux/scp/development/$BROAD_USER/mongo/hostname

export DATABASE_HOST="<ip from Vault (omit brackets)>"
```

If you are developing updates for Sentry logging, then set the DSN:

```
vault read secret/kdux/scp/production/scp_config.json | grep SENTRY

export SENTRY_DSN="<Sentry DSN value from Vault>"
```

Be sure to `unset SENTRY_DSN` when your updates are done, so development logs are not always sent to Sentry.

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

After [installing](#Install):

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

For more, see https://docs.pytest.org/en/stable/usage.html.

# Use

Run this every time you start a new terminal to work on this project:

```
source env/bin/activate
```

See [`ingest_pipeline.py`](https://github.com/broadinstitute/scp-ingest-pipeline/blob/ew-tests-hook/ingest/ingest_pipeline.py) for usage examples.
