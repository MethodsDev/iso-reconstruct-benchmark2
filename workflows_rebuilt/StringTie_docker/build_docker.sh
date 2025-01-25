#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:${VERSION} .
docker build -t us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:latest .


# verify installation
docker run --rm -it  us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:latest stringtie --version

