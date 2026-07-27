#!/usr/bin/env bash

set -euo pipefail
: "${SAI_SYSTEM_MPIRUN:?}"
: "${MAP_OPT:?}"

exec "$SAI_SYSTEM_MPIRUN" --map-by "$MAP_OPT" "$@"
