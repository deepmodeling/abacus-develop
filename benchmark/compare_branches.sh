#!/bin/bash
# Cross-branch PPCG benchmark comparison.
# Compares PPCG performance between two git branches.
#
# Usage: ./compare_branches.sh [base_branch] [target_branch] [--quick]
#   base_branch   — baseline branch (default: master)
#   target_branch — optimized branch (default: HEAD / current branch)
#   --quick       — use smaller matrix set

set -e

BASE_BRANCH="${1:-master}"
TARGET_BRANCH="${2:-HEAD}"
QUICK=""
if [[ "$3" == "--quick" ]] || [[ "$1" == "--quick" ]]; then
    QUICK="--quick"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MPIRUN=/opt/intel/oneapi/mpi/2021.13/bin/mpirun

echo "=== PPCG Cross-Branch Benchmark ==="
echo "Base:    $BASE_BRANCH"
echo "Target:  $TARGET_BRANCH"
echo ""

ORIG_BRANCH=$(cd "$REPO_DIR" && git branch --show-current)
STASHED=0

cleanup() {
    echo ""
    echo "=== Restoring original state ==="
    cd "$REPO_DIR"
    if git branch --show-current != "$ORIG_BRANCH" 2>/dev/null; then
        git checkout "$ORIG_BRANCH" 2>/dev/null || true
    fi
    if [ $STASHED -eq 1 ]; then
        git stash pop 2>/dev/null || true
    fi
}
trap cleanup EXIT

# Save any uncommitted changes
cd "$REPO_DIR"
if ! git diff-index --quiet HEAD -- 2>/dev/null; then
    git stash push -m "bench_compare_autostash" 2>/dev/null || true
    STASHED=1
fi

# Build and benchmark on base branch
echo "=== Benchmarking base branch: $BASE_BRANCH ==="
git checkout "$BASE_BRANCH" 2>/dev/null
CC=/opt/intel/oneapi/mpi/2021.13/bin/mpicc \
CXX=/opt/intel/oneapi/mpi/2021.13/bin/mpicxx \
cmake -B build -DBUILD_TESTING=ON -DENABLE_MPI=ON -DENABLE_LCAO=ON > /dev/null 2>&1
cmake --build build -j$(nproc) --target MODULE_HSOLVER_ppcg_bench > /dev/null 2>&1
bash "$SCRIPT_DIR/bench_ppcg.sh" $QUICK before.csv

# Build and benchmark on target branch
echo ""
echo "=== Benchmarking target branch: $TARGET_BRANCH ==="
git checkout "$TARGET_BRANCH" 2>/dev/null
CC=/opt/intel/oneapi/mpi/2021.13/bin/mpicc \
CXX=/opt/intel/oneapi/mpi/2021.13/bin/mpicxx \
cmake -B build -DBUILD_TESTING=ON -DENABLE_MPI=ON -DENABLE_LCAO=ON > /dev/null 2>&1
cmake --build build -j$(nproc) --target MODULE_HSOLVER_ppcg_bench > /dev/null 2>&1
bash "$SCRIPT_DIR/bench_ppcg.sh" $QUICK after.csv

# Generate comparison report
echo ""
echo "=== Comparison Report ==="
echo ""

if [ -f before.csv ] && [ -f after.csv ]; then
    echo "Configuration               Before(ms)  After(ms)  Speedup  Before(iter)  After(iter)"
    echo "---------------------------------------------------------------------------------------"

    # Skip header line of before.csv
    tail -n +2 before.csv | while IFS=, read -r npw nband sparsity mpi omp iter time err; do
        after_line=$(grep "^${npw},${nband},${sparsity},${mpi},${omp}," after.csv 2>/dev/null || echo "")
        if [ -n "$after_line" ]; then
            after_time=$(echo "$after_line" | cut -d, -f7)
            after_iter=$(echo "$after_line" | cut -d, -f6)
            if [ -n "$after_time" ] && [ -n "$time" ]; then
                speedup=$(echo "scale=2; $time / $after_time" | bc 2>/dev/null || echo "N/A")
                printf "%-28s %10.1f  %9.1f  %7s  %12s  %11s\n" \
                    "${npw}x${npw}/${nband}/s${sparsity}/mpi${mpi}/omp${omp}" \
                    "$time" "$after_time" "${speedup}x" "$iter" "$after_iter"
            fi
        fi
    done
    echo ""
    echo "Before results: before.csv"
    echo "After results:  after.csv"
else
    echo "Missing result files — benchmark may have failed."
fi
