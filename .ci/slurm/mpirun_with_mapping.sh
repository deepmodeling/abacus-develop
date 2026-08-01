#!/usr/bin/env bash

set -euo pipefail
: "${CI_SYSTEM_MPIRUN:?}"
: "${MAP_OPT:?}"

exec "$CI_SYSTEM_MPIRUN" --map-by "$MAP_OPT" "$@"
