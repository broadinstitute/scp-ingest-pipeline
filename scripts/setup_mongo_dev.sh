#! /bin/bash

# PREREQUISITE
# Set up your GitHub token in a file for use to access Vault
#
# To set up your development environment before running ingest pipeline on the command line, run: 
# source setup_mongo_dev.sh <path to your GitHub token file>
#
# Keep "Dev env vars" synced with `ingest-local-setup.bash`

VAULT_TOKEN_PATH="$1"
if [[ -z "$VAULT_TOKEN_PATH" ]]
then
  echo "You must provide a path to a GitHub token to proceed"
  exit 1
fi
vault login -method=github token=$(cat $VAULT_TOKEN_PATH)
if [[ $? -ne 0 ]]
then
  echo "Unable to authenticate into vault"
  exit 1
fi

# Dev env vars
export BROAD_USER=`whoami`
export MONGODB_USERNAME='single_cell'
export DATABASE_NAME='single_cell_portal_development'
export MONGODB_PASSWORD=`vault read secret/kdux/scp/development/$BROAD_USER/mongo/user | grep password | awk '{ print $2 }' `
export DATABASE_HOST=`vault read secret/kdux/scp/development/$BROAD_USER/mongo/hostname | grep ip | awk '{ $2=substr($2,2,length($2)-2); print $2 }' ` 
export BYPASS_MONGO_WRITES='yes'
export BARD_HOST_URL="https://terra-bard-dev.appspot.com"
