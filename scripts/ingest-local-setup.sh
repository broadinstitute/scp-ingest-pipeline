#! /bin/bash

# Called in docker-compose-setup.sh for Docker development environment
#
# Keep "Dev env vars" synced with `setup-mongo-dev.sh`

GOOGLE_PROJECT=$(gcloud info --format="value(config.project)")

# Dev env vars
BROAD_USER=`whoami`
MONGODB_USERNAME='single_cell'
DATABASE_NAME='single_cell_portal_development'
MONGODB_PASSWORD=$(gcloud secrets versions access latest --project=$GOOGLE_PROJECT --secret=mongo-user | jq .password)
DATABASE_HOST=$(gcloud secrets versions access latest --project=$GOOGLE_PROJECT --secret=mongo-hostname| jq -r '.ip[0]')
BYPASS_MONGO_WRITES='yes'
BARD_HOST_URL="https://terra-bard-dev.appspot.com"

mkdir -p config
echo "export BROAD_USER=${BROAD_USER}" >| config/.ingest_env
echo "export MONGODB_USERNAME=${MONGODB_USERNAME}" >> config/.ingest_env
echo "export DATABASE_NAME=${DATABASE_NAME}" >> config/.ingest_env
echo "export MONGODB_PASSWORD=${MONGODB_PASSWORD}" >> config/.ingest_env
echo "export DATABASE_HOST=${DATABASE_HOST}" >> config/.ingest_env
echo "export BYPASS_MONGO_WRITES=${BYPASS_MONGO_WRITES}" >> config/.ingest_env
echo "export BARD_HOST_URL=${BARD_HOST_URL}" >> config/.ingest_env
