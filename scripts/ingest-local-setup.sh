#! /bin/bash

## PREREQUISITE
## Set up your GitHub token in a file for use to access Vault

## To set up your development environment before running Ingest Pipeline on the command line, run: 
## source setup_mongo_dev.sh <path to your GitHub token file>

# eweitz 2023-01: Change these export to echo then pipe to file at path:
# config/secrets/.source_env.bash (do appropriate `mkdir -p` etc. before if needed)
# Add above path to .gitignore 

VAULT_TOKEN_PATH="$1"
if [[ -z "$VAULT_TOKEN_PATH" ]]
then
  echo "You must provide a path to a GitHub token to proceed"
  exit 1
fi
vault login -method=github token=$(cat $VAULT_TOKEN_PATH)
if [[ $? -ne 0 ]]
then
  echo "Unable to authenticate into Vault"
  exit 1
fi

BROAD_USER=`whoami`
MONGODB_USERNAME='single_cell'
DATABASE_NAME='single_cell_portal_development'
MONGODB_PASSWORD=`vault read secret/kdux/scp/development/$BROAD_USER/mongo/user | grep password | awk '{ print $2 }' `
DATABASE_HOST=`vault read secret/kdux/scp/development/$BROAD_USER/mongo/hostname | grep ip | awk '{ $2=substr($2,2,length($2)-2); print $2 }' ` 
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
