#!/bin/bash

# Based on CP2K's clang_format_wrapper.sh, adapted for ABACUS.

if (($# != 1)); then
  echo "Usage: clang_format_wrapper.sh <file>"
  exit 1
fi

STYLE_FILE="/opt/abacus-precommit/abacus_clang_format.yaml"
if [[ -f "${STYLE_FILE}" ]]; then
  CLANG_FORMAT_STYLE="file:${STYLE_FILE}"
else
  CLANG_FORMAT_STYLE="llvm"
fi

# For some files clang-format can require too much memory. As a work-around we
# use the "clang-format off" marker to disable formatting. To actually reduce
# the memory usage this script hides the offending sections from clang-format.
if csplit --quiet --prefix="section" "$1" "/clang-format off/"; then
  clang-format --style="${CLANG_FORMAT_STYLE}" -i section00
  cat section* > "$1"
  rm section*
else
  clang-format --style="${CLANG_FORMAT_STYLE}" -i "$1"
fi

# EOF
