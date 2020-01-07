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
FROM marketplace.gcr.io/google/ubuntu1804:latest

# RUN echo "Uncomment to clear cached layers below this statement (2020-01-07-0947)"

# Install Python 3.7
RUN apt-get -y update && \
  apt -y install software-properties-common && \
  add-apt-repository ppa:deadsnakes/ppa && \
  apt -y install python3-pip && \
  apt -y install python3.7

RUN python3.7 -m pip install pip

# Set cleaner defaults (`alias` fails)
RUN ln -s /usr/bin/python3.7 /usr/bin/python && \
  ln -s /usr/bin/pip3 /usr/bin/pip

# Copy contents of this repo into the Docker image
# (See .Dockerignore for omitted files)
COPY . scp-ingest-pipeline

WORKDIR /scp-ingest-pipeline

# Install Python dependencies
RUN python3.7 -m pip install -r requirements.txt

WORKDIR /scp-ingest-pipeline/ingest
CMD ["python", "ingest_pipeline.py", "--help"]