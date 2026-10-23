#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 VALIDATION_TOOL_DIRECTORY" >&2
    exit 2
fi

python3 -m unittest discover -s "$1" -p 'test_librpa_*.py'
