# Dockerfile for SCP Ingest Service

# Use a managed base image from Google.  It is continually updated for
# security fixes (thus the "latest" tag).
# https://github.com/GoogleContainerTools/base-images-docker/tree/master/ubuntu
FROM marketplace.gcr.io/google/ubuntu1804:latest

# Install Python 3.7
RUN apt update
RUN apt -y install software-properties-common dirmngr apt-transport-https lsb-release ca-certificates
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt -y install python3.7
RUN alias python=python3.7

# Copy contents of this repo into the Docker image
# (See .Dockerignore for omitted files)
ADD . scp-ingest-service