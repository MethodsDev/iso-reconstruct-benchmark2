#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest .

# verify installation
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:${VERSION} IsoQuant-runner.py -h


