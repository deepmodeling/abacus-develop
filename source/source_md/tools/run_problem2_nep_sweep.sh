#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MD_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DIR="$(cd "${MD_DIR}/../.." && pwd)"

CXX="${CXX:-}"
if [[ -z "${CXX}" ]]; then
  if command -v g++-15 >/dev/null 2>&1; then
    CXX="g++-15"
  elif command -v g++ >/dev/null 2>&1; then
    CXX="g++"
  else
    echo "No C++ compiler found. Set CXX to a compiler with OpenMP support." >&2
    exit 1
  fi
fi

MODEL_FILE="${1:-${REPO_DIR}/tests/PP_ORB/nep_hfo2.txt}"
SIZES="${2:-64,216,512}"
THREADS="${3:-1,2,4}"
REPEATS="${4:-10}"
OUT_CSV="${5:-${MD_DIR}/reports/problem2_nep_sweep_results.csv}"
BUILD_DIR="${MD_DIR}/build_problem2_nep_sweep"
BIN="${BUILD_DIR}/nep_omp_sweep"

mkdir -p "${BUILD_DIR}" "$(dirname "${OUT_CSV}")"

"${CXX}" -std=c++14 -O2 -fopenmp \
  -I"${MD_DIR}" \
  "${MD_DIR}/tools/nep_omp_sweep.cpp" \
  "${MD_DIR}/potential/ml/nep/nep_cpu.cpp" \
  "${MD_DIR}/potential/ml/nep/neighbor_nep.cpp" \
  "${MD_DIR}/potential/ml/nep/ewald_nep.cpp" \
  -o "${BIN}"

"${BIN}" "${MODEL_FILE}" "${SIZES}" "${THREADS}" "${REPEATS}" | tee "${OUT_CSV}"

echo "Saved CSV to ${OUT_CSV}" >&2
