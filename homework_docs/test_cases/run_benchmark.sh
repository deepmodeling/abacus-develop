#!/bin/bash
# Workflow C Baseline Benchmark Script
# Runs 3 test cases x 6 MPI/OMP configs x 3 repetitions = 54 runs
set -e

ABACUS_BIN="/home/aunixt/abacus-develop/build/abacus_basic_para"
TEST_BASE="/home/aunixt/abacus-develop/homework_docs/test_cases"
RESULT_BASE="/home/aunixt/abacus-develop/homework_docs/test_cases/results"

mkdir -p "$RESULT_BASE"

# Test configurations: case_name directory
CASES=("gaas_tiny" "gaas_small" "gaas_medium")

# MPI x OMP matrix
# Format: "np omp label"
CONFIGS=(
    "1 1 serial"
    "1 4 omp4"
    "1 8 omp8"
    "2 2 mix_np2_omp2"
    "2 4 mix_np2_omp4"
    "4 2 mix_np4_omp2"
)

REPS=3
TIMEOUT_PER_RUN=600  # 10 minutes per run

for case_name in "${CASES[@]}"; do
    CASE_DIR="$TEST_BASE/$case_name"
    CASE_RESULT="$RESULT_BASE/$case_name"
    mkdir -p "$CASE_RESULT"

    echo "============================================="
    echo "Testing: $case_name"
    echo "============================================="

    for config in "${CONFIGS[@]}"; do
        read -r np omp label <<< "$config"

        for rep in $(seq 1 $REPS); do
            run_name="np${np}_omp${omp}_r${rep}"
            stdout_file="$CASE_RESULT/${run_name}.stdout"
            stderr_file="$CASE_RESULT/${run_name}.stderr"

            echo "  [$case_name] $label (np=$np omp=$omp) rep=$rep ..."

            cd "$CASE_DIR"
            export OMP_NUM_THREADS=$omp

            if [ "$np" -eq 1 ]; then
                timeout $TIMEOUT_PER_RUN "$ABACUS_BIN" > "$stdout_file" 2> "$stderr_file" || {
                    echo "    WARNING: Run $run_name exited with code $?" | tee -a "$CASE_RESULT/errors.log"
                }
            else
                timeout $TIMEOUT_PER_RUN mpirun --allow-run-as-root -np $np "$ABACUS_BIN" > "$stdout_file" 2> "$stderr_file" || {
                    echo "    WARNING: Run $run_name exited with code $?" | tee -a "$CASE_RESULT/errors.log"
                }
            fi

            # Quick check for timer output
            if grep -q "TIME STATISTICS" "$stdout_file"; then
                echo "    OK (timer data found)"
            else
                echo "    WARNING: No TIME STATISTICS in $run_name" | tee -a "$CASE_RESULT/errors.log"
            fi

            # Clean up output directory to save space
            rm -rf OUT.ABACUS
        done
    done
done

echo ""
echo "============================================="
echo "Benchmark complete!"
echo "Results in: $RESULT_BASE"
echo "============================================="

# Quick summary
for case_name in "${CASES[@]}"; do
    CASE_RESULT="$RESULT_BASE/$case_name"
    total_runs=$(ls "$CASE_RESULT"/*.stdout 2>/dev/null | wc -l)
    timer_runs=$(grep -l "TIME STATISTICS" "$CASE_RESULT"/*.stdout 2>/dev/null | wc -l)
    echo "  $case_name: $timer_runs/$total_runs runs with timer data"
done
