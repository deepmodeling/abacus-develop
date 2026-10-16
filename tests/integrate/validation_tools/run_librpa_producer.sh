#!/bin/bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 ABACUS_EXECUTABLE CASES_DIRECTORY" >&2
    exit 2
fi

abacus=$1
cases_dir=$2
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

cd "$cases_dir"
bash "$script_dir/../Autotest.sh" -a "$abacus" -n 1 -j auto -p -f CASES_LIBRPA_PRODUCER.txt

while read -r case_name; do
    [ -z "$case_name" ] && continue
    case_dir="$cases_dir/$case_name"
    (cd "$case_dir" && python3 "$script_dir/check_librpa_producer.py" \
        --root . --manifest librpa_producer_manifest.json)
done < CASES_LIBRPA_PRODUCER.txt
