#! /bin/sh

# docker-compose-setup.sh
# bring up local development environment via docker-compose
# More context: 



usage=$(
cat <<EOF
$0 [OPTION]
-d   run docker-compose in detached mode (default is attatched to terminal STDOUT)
-c   enable VITE_FRONTEND_SERVICE_WORKER_CACHE (default is disabled)
-h   print this text
EOF
)

VAULT_TOKEN_PATH=""
while getopts "i:t:h" OPTION; do
case $OPTION in
  i)
    echo "### SETTING GCR IMAGE ###"
    export GCR_IMAGE="$OPTARG"
    ;;
  t)
    echo "### SETTING VAULT TOKEN ###"
    VAULT_TOKEN_PATH="$OPTARG"
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

if [[ $GCR_IMAGE = "" ]]; then
  IMAGE_NAME="gcr.io/broad-singlecellportal-staging/scp-ingest-pipeline-development"
  LATEST_TAG=$(gcloud container images list-tags ${IMAGE_NAME} --format='get(tags)' | head -n 1)
  export GCR_IMAGE="${IMAGE_NAME}:${LATEST_TAG}"
fi


if [[ $VAULT_TOKEN_PATH = "" ]]; then
  echo "Did not provide VAULT_TOKEN_PATH"
  exit 1
fi

echo "### SETTING UP ENVIRONMENT ###"
./scripts/ingest_local_setup.bash $VAULT_TOKEN_PATH

# If updating Python dependencies via pip
# comment out `docker pull` line below, commit 
docker pull $GCR_IMAGE

# Run the "ingest" service (define in docker-compose-dev.yaml), 
# then enter the Bash shell
docker-compose -f docker-compose-dev.yaml run ingest bash
