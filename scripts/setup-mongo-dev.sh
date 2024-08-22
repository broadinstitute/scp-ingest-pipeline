#! /bin/bash

# To set up your development environment before running ingest pipeline on the command line, run:
# source setup-mongo-dev.sh
#
# Keep "Dev env vars" synced with `ingest-local-setup.bash`

GOOGLE_PROJECT=$(gcloud info --format="value(config.project)")
MONGODB_PASSWORD=$(gcloud secrets versions access latest --project=$GOOGLE_PROJECT --secret=mongo-user | jq .password)
DATABASE_HOST=$(gcloud secrets versions access latest --project=$GOOGLE_PROJECT --secret=mongo-hostname | jq -r '.ip[0]')

# Dev env vars
export MONGODB_USERNAME='single_cell'
export DATABASE_NAME='single_cell_portal_development'
export MONGODB_PASSWORD=$MONGODB_PASSWORD
export DATABASE_HOST=$DATABASE_HOST
export BYPASS_MONGO_WRITES='yes'
export BARD_HOST_URL="https://terra-bard-dev.appspot.com"
