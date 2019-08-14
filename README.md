# scp-ingest-pipeline
File Ingest Pipeline for Single Cell Portal

[![Build status](https://img.shields.io/circleci/build/github/broadinstitute/scp-ingest-pipeline.svg)](https://circleci.com/gh/broadinstitute/scp-ingest-pipeline)

The SCP Ingest Pipeline is an ETL pipeline for single-cell RNA-seq data.  

# Prerequisites
* Python 3.6+
* Google Cloud project, including
  * Firestore in Native Mode
  * Service accounts

# Install
```
git clone https://github.com/broadinstitute/scp-ingest-pipeline
python3 -m venv env --copies
pip install -r requirements.txt
```

# Test
After installing:
```
cd tests
pytest
```

# Commit and push
Install pre-push hooks to automatically ensure that Ingest Pipeline tests pass:
```
pre-commit install -t pre-push
```