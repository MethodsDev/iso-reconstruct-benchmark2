#!/bin/bash
set -ex

VERSION=$(cat VERSION.txt)
IMAGE=us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2

docker push ${IMAGE}:${VERSION}
docker push ${IMAGE}:latest
