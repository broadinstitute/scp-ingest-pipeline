# Dockerfile for SCP Ingest Service
#
# PREREQUISITES
# Run the following command to register gcloud as a Docker credential helper:
# gcloud auth configure-docker

# TODO:
# ^ Put this in a Bash script

# Use a managed base image from Google.  It is continually updated for
# security fixes (thus the "latest" tag).
# https://github.com/GoogleContainerTools/base-images-docker/tree/master/ubuntu
FROM marketplace.gcr.io/google/ubuntu1804:latest

RUN echo "Uncomment to clear cached layers below this statement (20190703-1646)"

# Install Python 3.6
# (Auxiliary Ubuntu packages for Python assume Python 3.6;
# configuring them for Python >= 3.7 is not straightforward.)
RUN apt-get -y update
RUN apt -y install software-properties-common dirmngr apt-transport-https lsb-release ca-certificates
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt -y install python3.6
RUN apt -y install python3-pip

# Set cleaner defaults (`alias` fails)
# These `ln` commands do not persist when running Bash in container -- why?
RUN ln -s /usr/bin/python3 /usr/bin/python & \
    ln -s /usr/bin/pip3 /usr/bin/pip

# Copy contents of this repo into the Docker image
# (See .Dockerignore for omitted files)
COPY . scp-ingest-service

WORKDIR /scp-ingest-service

RUN pip install -r requirements.txt

WORKDIR /scp-ingest-service/ingest
CMD ["python3", "ingest.py", "--help"]