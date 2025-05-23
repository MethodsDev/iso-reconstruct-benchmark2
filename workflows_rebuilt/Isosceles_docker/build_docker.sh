#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:latest .

# verify installation
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:${VERSION}

