#!/bin/bash

# =========================================================================
# MPI Parallel Optimization Test Script
# =========================================================================
# This script runs the MPI unit tests for the eigenvalue solver with
# different numbers of processes to verify:
#   - Correctness across process counts
#   - Performance scaling
#   - Communication error handling
#
# Usage: ./diago_mpi_parallel_test.sh
# =========================================================================

set -e

# Detect number of available cores
np=$(cat /proc/cpuinfo 2>/dev/null | grep "cpu cores" | uniq | awk '{print $NF}' || echo 4)
echo "[INFO] Available cores: $np"

# Test executable name
TEST_EXE="./MODULE_HSOLVER_mpi"
if [ ! -f "$TEST_EXE" ]; then
    echo "[ERROR] Test executable $TEST_EXE not found"
    echo "[INFO] Please build with: cmake --build . --target MODULE_HSOLVER_mpi"
    exit 1
fi

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track results
PASS_COUNT=0
FAIL_COUNT=0
TOTAL_TESTS=0

# =========================================================================
# Function: run_mpi_test
# =========================================================================
run_mpi_test() {
    local nprocs=$1
    local label=$2

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    echo ""
    echo "============================================================"
    echo "[TEST] $label (nprocs=$nprocs)"
    echo "============================================================"

    if OMP_NUM_THREADS=1 mpirun --allow-run-as-root -np "$nprocs" "$TEST_EXE" 2>&1; then
        echo -e "${GREEN}[  PASSED  ]${NC} $label with $nprocs processes"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        echo -e "${RED}[  FAILED  ]${NC} $label with $nprocs processes"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
}

# =========================================================================
# Test with different process counts
# =========================================================================

echo "============================================================"
echo " MPI Parallel Eigenvalue Solver Optimization Test Suite"
echo "============================================================"
echo ""

# Determine which process counts to test
# Test at least 1, 2, 3, 4 (or min(nprocs, 1..4))

for nproc in 1 2 3 4; do
    if [ "$nproc" -le "$np" ]; then
        run_mpi_test "$nproc" "MPI Correctness ($nproc procs)"
    fi
done

# Additional test with more processes if available
if [ "$np" -ge 6 ]; then
    run_mpi_test 6 "MPI Correctness (6 procs)"
fi

if [ "$np" -ge 8 ]; then
    run_mpi_test 8 "MPI Correctness (8 procs)"
fi

# =========================================================================
# Summary
# =========================================================================

echo ""
echo "============================================================"
echo " Test Summary"
echo "============================================================"
echo -e "Total:  $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASS_COUNT${NC}"
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo -e "${RED}Failed: $FAIL_COUNT${NC}"
else
    echo -e "Failed: $FAIL_COUNT"
fi
echo "============================================================"

if [ "$FAIL_COUNT" -gt 0 ]; then
    echo -e "${RED}[FAIL] Some MPI tests failed!${NC}"
    exit 1
else
    echo -e "${GREEN}[PASS] All MPI tests passed!${NC}"
    exit 0
fi
