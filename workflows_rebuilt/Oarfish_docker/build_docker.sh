#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish:latest .

# verify installation
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish Oarfish_runner.py -h


