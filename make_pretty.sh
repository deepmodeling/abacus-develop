#!/bin/bash -e

ABACUS_ROOT=$(dirname "$0")

"${ABACUS_ROOT}/tools/precommit/precommit.py" --allow-modifications "$@"

#EOF
