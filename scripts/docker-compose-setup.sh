#! /bin/sh

# docker-compose-setup.sh
# bring up local development environment via docker-compose
#
# More context:
# See https://github.com/broadinstitute/scp-ingest-pipeline#docker

usage=$(
cat <<EOF
$0 [OPTION]
-i   Set URL for GCR image; helpful if not using latest development
-h   print this text
EOF
)

GCR_IMAGE=""
VAULT_TOKEN_PATH=""
while getopts "i:h" OPTION; do
case $OPTION in
  i)
    echo "### SETTING GCR IMAGE ###"
    export GCR_IMAGE="$OPTARG"
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

echo "### SETTING UP ENVIRONMENT ###"
./scripts/ingest-local-setup.sh

docker pull $GCR_IMAGE

# Run the "ingest" service (define in docker-compose-dev.yaml),
# then enter the Bash shell
docker-compose -f docker-compose-dev.yaml run ingest bash
