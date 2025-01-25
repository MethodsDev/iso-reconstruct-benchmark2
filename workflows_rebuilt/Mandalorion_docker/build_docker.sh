#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion:latest .


# verify installation
docker run us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion:${VERSION} Mandalorian-runner.py -h

