#!/bin/bash
set -ex

VERSION=$(cat VERSION.txt)

IMAGE=us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2

docker build -t ${IMAGE}:${VERSION} .
docker build -t ${IMAGE}:latest .

docker run --rm ${IMAGE}:latest --help
