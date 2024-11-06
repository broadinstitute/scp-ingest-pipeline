#! /bin/bash

# vault-gsm-migration.sh
# script to pull out a single ingest-related secret from vault and import it into Google Secret Manager (GSM)
# requires gcloud, vault, and jq to be installed and configured
# adapted from https://github.com/broadinstitute/single_cell_portal_core/blob/development/bin/vault-gsm-migration.sh

function authenticate_vault {
  TOKEN_PATH="$1"
  vault login -method=github token="$(cat "$TOKEN_PATH")"
}

function extract_secret_from_vault {
  SECRET_PATH="$1"
  vault read -field=data -format=json "$SECRET_PATH" | jq
}

function delete_gsm_secret {
  SECRET_NAME="$1"
  gcloud secrets delete "$SECRET_NAME" --quiet
}

function copy_secret_to_gsm {
  VAULT_PATH="$1"
  SECRET_NAME="$2"
  GOOGLE_PROJECT="$3"
  extract_secret_from_vault "$VAULT_PATH" | gcloud secrets create "$SECRET_NAME" --project="$GOOGLE_PROJECT" \
                                            --replication-policy=automatic --data-file=-
}

function exit_with_error {
  echo "### ERROR: $1 ###"
  exit 1
}

function main {
  GOOGLE_PROJECT=$(gcloud info --format="value(config.project)")
  BASEPATH="secret/kdux/scp/development/$(whoami)"
  VAULT_SECRET=""
  ROLLBACK="false"
  NAME_OVERRIDE=""
  while getopts "t:p:b:s:n:rh" OPTION; do
    case $OPTION in
      t)
        VAULT_TOKEN="$OPTARG"
        ;;
      p)
        GOOGLE_PROJECT="$OPTARG"
        echo "setting GOOGLE_PROJECT to $GOOGLE_PROJECT"
        ;;
      b)
        BASEPATH="$OPTARG"
        echo "setting BASEPATH to $BASEPATH"
        ;;
      s)
        VAULT_SECRET="$OPTARG"
        echo "setting VAULT_SECRET to $VAULT_SECRET"
        ;;
      n)
        NAME_OVERRIDE="$OPTARG"
        echo "using $NAME_OVERRIDE as GSM secret name"
        ;;
      r)
        ROLLBACK="true"
        ;;
      h)
        echo "$usage"
        exit 0
        ;;
      *)
        echo "unrecognized option"
        echo "$usage"
        exit 1
        ;;
    esac
  done

  if [[ -z "$VAULT_SECRET" ]]; then
    echo "ERROR: Did not provide relative path to vault secret to copy"
    echo "$usage"
    exit 1
  fi

  # set up for GSM secret based on existing vault path
  # changes / to -, e.g. mongo/hostname becomes mongo-hostname
  GSM_SECRET_NAME=$(echo "$VAULT_SECRET" | sed 's/[/_]/\-/g')

  # use name override if specified
  if [[ -n "$NAME_OVERRIDE" ]]; then
    GSM_SECRET_NAME="$NAME_OVERRIDE"
  fi

  if [[ "$ROLLBACK" = "true" ]]; then
    echo "Rolling back GSM migration"
    echo -n "Deleting $VAULT_SECRET... "
    delete_gsm_secret "$GSM_SECRET_NAME"
    echo "Rollback complete!"
    exit 0
  fi

  # construct vault path
  VAULT_SECRET_PATH="$BASEPATH/$VAULT_SECRET"

  authenticate_vault "$VAULT_TOKEN" || exit_with_error "cannot authenticate into vault"

  echo -n "Copying $VAULT_SECRET from vault... "
  copy_secret_to_gsm "$VAULT_SECRET_PATH" "$GSM_SECRET_NAME" "$GOOGLE_PROJECT" || exit_with_error "failed to copy $VAULT_SECRET"

  echo "Secret migrated to GSM"
  exit 0
}

usage=$(
cat <<EOF
USAGE:
  $(basename "$0") -t VAULT_TOKEN -s VAULT_SECRET [<options>]

  -t VAULT_TOKEN  set the path to token used to authenticate into vault (no default)
  -s VAULT_SECRET set the name of the requested secret in vault to copy (no default)

  [OPTIONS]
  -n GSM_NAME     set the name of the new GSM secret, otherwise defaults to old name with special chars converted to -
  -p PROJECT      set the GCP project in which to create secrets (defaults to '$(gcloud info --format="value(config.project)")')
  -b BASEPATH     set the base vault path for retrieving secrets (defaults to 'secret/kdux/scp/development/$(whoami)')
  -r              roll back migration and delete scp-ingest-pipeline secret in GSM
  -h              print this message
EOF
)

main "$@"
