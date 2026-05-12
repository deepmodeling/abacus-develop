#!/bin/bash -e

# Based on CP2K's tools/precommit/start_local_server.sh, adapted for ABACUS.

set -x

SHORT_SHA=$(git rev-parse --short HEAD)
podman build --build-arg "REVISION=${SHORT_SHA}" -t abacus-precommit .
podman run --rm -p127.0.0.1:8080:8080 abacus-precommit

# EOF
