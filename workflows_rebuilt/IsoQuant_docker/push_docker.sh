#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

docker push us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:${VERSION}
docker push us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest 


