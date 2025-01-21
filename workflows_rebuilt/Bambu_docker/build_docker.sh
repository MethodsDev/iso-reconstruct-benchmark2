#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest .

# verify installation
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:${VERSION} bambu-runner.Rscript

