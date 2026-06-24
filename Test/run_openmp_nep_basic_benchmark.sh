#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
RESULT_DIR="${SCRIPT_DIR}/results"
BIN="${BUILD_DIR}/openmp_nep_basic_benchmark"
CSV="${RESULT_DIR}/openmp_nep_basic_benchmark.csv"
LOG="${RESULT_DIR}/openmp_nep_basic_benchmark.log"

NATOMS="${NATOMS:-2000000}"
REPEAT="${REPEAT:-5}"
CXX="${CXX:-g++}"

mkdir -p "${BUILD_DIR}" "${RESULT_DIR}"

{
  echo "Compiler: $(${CXX} --version | head -n 1)"
  echo "NATOMS=${NATOMS}"
  echo "REPEAT=${REPEAT}"
  echo "Build: ${CXX} -O3 -std=c++17 -fopenmp"
} > "${LOG}"

"${CXX}" -O3 -std=c++17 -fopenmp "${SCRIPT_DIR}/openmp_nep_basic_benchmark.cpp" -o "${BIN}" 2>&1 | tee -a "${LOG}"

: > "${CSV}"
for threads in 1 2 4 8; do
  echo "Running with OMP_NUM_THREADS=${threads}" | tee -a "${LOG}"
  export OMP_NUM_THREADS="${threads}"
  export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
  export OMP_PLACES="${OMP_PLACES:-cores}"
  tmp_csv="${RESULT_DIR}/run_${threads}.csv"
  "${BIN}" --threads "${threads}" --natoms "${NATOMS}" --repeat "${REPEAT}" > "${tmp_csv}"
  if [[ "${threads}" == "1" ]]; then
    cat "${tmp_csv}" >> "${CSV}"
  else
    tail -n +2 "${tmp_csv}" >> "${CSV}"
  fi
done

echo "CSV: ${CSV}" | tee -a "${LOG}"
