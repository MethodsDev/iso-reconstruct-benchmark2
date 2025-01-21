#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`


docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:${VERSION} .
docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest .

