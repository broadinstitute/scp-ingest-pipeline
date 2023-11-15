# Dockerfile for SCP Ingest Pipeline
#
# PREREQUISITES
# Run the following command to register gcloud as a Docker credential helper:
# gcloud auth configure-docker

# TODO:
# ^ Put this in a Bash script

# Use a managed base image from Google.  It is continually updated for
# security fixes (thus the "latest" tag).
# https://github.com/GoogleContainerTools/base-images-docker/tree/master/ubuntu
FROM marketplace.gcr.io/google/ubuntu2004:latest

# RUN echo "Uncomment to clear cached layers below this statement (2022-03-14-1441)"

# Install Python 3.10
RUN apt-get -y update && \
  apt-get -y install software-properties-common && \
  add-apt-repository ppa:deadsnakes/ppa && \
  apt-get -y install python3-pip && \
  apt-get -y install python3.10 && \
  apt-get -y install python3.10-dev && \
  apt-get -y install python3.10-distutils

RUN apt-get update && apt-get install curl

RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

# symlink python3.10 to python
RUN ln -s /usr/bin/python3.10 /usr/bin/python

# Copy contents of this repo into the Docker image
# (See .Dockerignore for omitted files)
COPY . scp-ingest-pipeline

WORKDIR /scp-ingest-pipeline

# Install Python dependencies
RUN python -m pip install -r requirements.txt

WORKDIR /scp-ingest-pipeline/ingest
CMD ["python", "ingest_pipeline.py", "--help"]
