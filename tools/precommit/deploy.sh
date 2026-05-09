#!/bin/bash -e

# Based on CP2K's tools/precommit/deploy.sh, adapted for ABACUS.
# Adjust PROJECT/REGISTRY/REGION before using for production deployment.

set -x

: "${ABACUS_PRECOMMIT_IMAGE:=abacus-precommit}"
: "${ABACUS_PRECOMMIT_REGISTRY:=gcr.io/abacus-project/img_abacus_precommit}"
: "${ABACUS_PRECOMMIT_REGION:=us-central1}"

SHORT_SHA=$(git rev-parse --short HEAD)
docker build --build-arg "REVISION=${SHORT_SHA}" -t "${ABACUS_PRECOMMIT_IMAGE}" .

docker tag "${ABACUS_PRECOMMIT_IMAGE}" "${ABACUS_PRECOMMIT_REGISTRY}:${SHORT_SHA}"
docker tag "${ABACUS_PRECOMMIT_IMAGE}" "${ABACUS_PRECOMMIT_REGISTRY}:latest"
docker push "${ABACUS_PRECOMMIT_REGISTRY}:${SHORT_SHA}"
docker push "${ABACUS_PRECOMMIT_REGISTRY}:latest"

gcloud run deploy abacus-precommit --platform=managed --region="${ABACUS_PRECOMMIT_REGION}" \
  --image="${ABACUS_PRECOMMIT_REGISTRY}:${SHORT_SHA}"

# EOF
