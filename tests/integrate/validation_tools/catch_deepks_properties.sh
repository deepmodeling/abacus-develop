#!/bin/bash

# Standalone entry for the DeePKS property collector.
#
# catch_properties.sh sources props_deepks.sh directly; this thin wrapper is
# kept for manual invocation, e.g. `bash catch_deepks_properties.sh result`.
# props_init() is intentionally NOT called: standalone mode must append to
# the caller-supplied result file without truncating it, and
# run_deepks_props() reads the few INPUT switches it needs itself.
PROPS_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$PROPS_SCRIPT_DIR/props_common.sh"
props_result_file=$1
source "$PROPS_SCRIPT_DIR/props_deepks.sh"
run_deepks_props
