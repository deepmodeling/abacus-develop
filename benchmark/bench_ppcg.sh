#!/bin/bash
# PPCG benchmark — measures runtime and iterations across matrix sizes and thread counts.
#
# Usage: ./bench_ppcg.sh [--quick] [output.csv]
#   --quick: smaller matrix set for fast validation

set -e

MPIRUN=/opt/intel/oneapi/mpi/2021.13/bin/mpirun
BUILD_DIR=$(cd "$(dirname "$0")/../build" && pwd)
BENCH_BIN="$BUILD_DIR/source/source_hsolver/test/MODULE_HSOLVER_ppcg_bench"

OUTPUT="${1:-ppcg_bench_results.csv}"

# Test configurations: npw nband sparsity ethr
if [[ "$1" == "--quick" ]]; then
    shift
    OUTPUT="${1:-ppcg_bench_results.csv}"
    CONFIGS=(
        "100  10  0  1e-7"
        "200  20  6  1e-7"
    )
else
    CONFIGS=(
        "100  10  0  1e-7"
        "500  50  6  1e-7"
        "1000 100 8  1e-7"
        "200  20  5  1e-7"  # closely spaced eigenvalues
    )
fi

OMP_THREADS=(1 2 4)

# CSV header (to stdout)
echo "npw,nband,sparsity,mpi_procs,omp_threads,iterations,time_ms,max_error"

for cfg in "${CONFIGS[@]}"; do
    read -r npw nband sparsity ethr <<< "$cfg"
    for omp in "${OMP_THREADS[@]}"; do
        export OMP_NUM_THREADS=$omp
        $MPIRUN -np 1 $BENCH_BIN $npw $nband $sparsity $ethr 2>/dev/null || echo "${npw},${nband},${sparsity},1,${omp},FAIL,FAIL,FAIL"
    done
done
