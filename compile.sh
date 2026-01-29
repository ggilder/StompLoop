#!/bin/bash

set -euo pipefail

docker build -t pd-armhf-build .
docker run --rm -it \
  --platform linux/arm/v7 \
  -v "$PWD:/work" \
  pd-armhf-build \
  bash -c "cd /work && make clean all"
