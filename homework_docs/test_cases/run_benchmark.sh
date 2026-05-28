#!/usr/bin/env bash
# Workflow C cache-reuse validation runner.
#
# Default behavior:
#   1. build the ABACUS binary and module_pw tests
#   2. run cache-related unit tests
#   3. run the selected SCF benchmark matrix
#   4. parse TIME STATISTICS into summary.txt
#
# Useful quick run:
#   ./run_benchmark.sh --quick
#
# Full run:
#   ./run_benchmark.sh

set -u
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
BUILD_DIR="$REPO_ROOT/build"
ABACUS_BIN="$BUILD_DIR/abacus_basic_para"
RESULT_ROOT="$SCRIPT_DIR/results_cache_$(date +"%Y%m%d_%H%M%S")"

CASES_CSV="gaas_tiny,gaas_small,gaas_medium"
CONFIGS_CSV="1:1:serial,1:4:omp4,1:8:omp8,2:2:mix_np2_omp2,2:4:mix_np2_omp4,4:2:mix_np4_omp2"
REPS=3
TIMEOUT_PER_RUN=600
DO_BUILD=1
DO_UNIT=1
DO_BENCH=1
DO_PARSE=1
ALLOW_FAILURES=1

usage() {
    cat <<'EOF'
Usage: run_benchmark.sh [options]

Options:
  --quick                 Run a short validation: gaas_tiny, serial + np2_omp2, 1 rep.
  --cases A,B,C           Case directories under homework_docs/test_cases.
  --configs spec          Comma list of np:omp:label entries.
  --reps N                Repetitions per case/config.
  --timeout SEC           Timeout per SCF run.
  --bin PATH              ABACUS executable.
  --build-dir PATH        CMake build directory.
  --out PATH              Result directory.
  --no-build              Skip cmake builds.
  --no-unit               Skip unit tests.
  --no-bench              Skip SCF benchmark runs.
  --no-parse              Skip parse_timers.py summary.
  --strict                Stop on first failed run.
  -h, --help              Show this help.

Every command is timestamped in run_log.tsv with start/end ISO time,
epoch-ns values, status, and elapsed seconds.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --quick)
            CASES_CSV="gaas_tiny"
            CONFIGS_CSV="1:1:serial,2:2:mix_np2_omp2"
            REPS=1
            ;;
        --cases)
            CASES_CSV="$2"
            shift
            ;;
        --configs)
            CONFIGS_CSV="$2"
            shift
            ;;
        --reps)
            REPS="$2"
            shift
            ;;
        --timeout)
            TIMEOUT_PER_RUN="$2"
            shift
            ;;
        --bin)
            ABACUS_BIN="$2"
            shift
            ;;
        --build-dir)
            BUILD_DIR="$2"
            ABACUS_BIN="$BUILD_DIR/abacus_basic_para"
            shift
            ;;
        --out)
            RESULT_ROOT="$2"
            shift
            ;;
        --no-build)
            DO_BUILD=0
            ;;
        --no-unit)
            DO_UNIT=0
            ;;
        --no-bench)
            DO_BENCH=0
            ;;
        --no-parse)
            DO_PARSE=0
            ;;
        --strict)
            ALLOW_FAILURES=0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

mkdir -p "$RESULT_ROOT"
LOG_FILE="$RESULT_ROOT/run_log.tsv"
SUMMARY_FILE="$RESULT_ROOT/summary.txt"
ERROR_FILE="$RESULT_ROOT/errors.log"

printf "phase\tcase\tconfig\trep\tcommand\tstart_iso\tend_iso\tstart_ns\tend_ns\telapsed_s\tstatus\n" > "$LOG_FILE"

now_iso() {
    date +"%Y-%m-%dT%H:%M:%S%z"
}

now_ns() {
    date +"%s%N"
}

elapsed_seconds() {
    awk -v start="$1" -v end="$2" 'BEGIN { printf "%.3f", (end - start) / 1000000000.0 }'
}

record_command() {
    local phase="$1"
    local case_name="$2"
    local config_label="$3"
    local rep="$4"
    local stdout_file="$5"
    local stderr_file="$6"
    shift 6

    local start_iso end_iso start_ns end_ns elapsed status command_text
    command_text="$*"
    start_iso=$(now_iso)
    start_ns=$(now_ns)

    echo "[$start_iso] START $phase case=$case_name config=$config_label rep=$rep"
    echo "  command: $command_text"

    "$@" > "$stdout_file" 2> "$stderr_file"
    status=$?

    end_ns=$(now_ns)
    end_iso=$(now_iso)
    elapsed=$(elapsed_seconds "$start_ns" "$end_ns")
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$phase" "$case_name" "$config_label" "$rep" "$command_text" \
        "$start_iso" "$end_iso" "$start_ns" "$end_ns" "$elapsed" "$status" >> "$LOG_FILE"

    echo "[$end_iso] END   $phase case=$case_name config=$config_label rep=$rep status=$status elapsed=${elapsed}s"

    if [[ "$status" -ne 0 ]]; then
        echo "FAILED phase=$phase case=$case_name config=$config_label rep=$rep status=$status" >> "$ERROR_FILE"
        if [[ "$ALLOW_FAILURES" -eq 0 ]]; then
            exit "$status"
        fi
    fi

    return "$status"
}

run_builds() {
    local out_dir="$RESULT_ROOT/build"
    mkdir -p "$out_dir"

    record_command build repo abacus 0 "$out_dir/abacus.stdout" "$out_dir/abacus.stderr" \
        cmake --build "$BUILD_DIR" --target abacus_basic_para -j2
    record_command build repo module_pw_serial 0 "$out_dir/module_pw_serial.stdout" "$out_dir/module_pw_serial.stderr" \
        cmake --build "$BUILD_DIR" --target MODULE_PW_basis_pw_k_serial -j2
    record_command build repo module_pw_mpi 0 "$out_dir/module_pw_mpi.stdout" "$out_dir/module_pw_mpi.stderr" \
        cmake --build "$BUILD_DIR" --target MODULE_PW_pw_test -j2
}

run_unit_tests() {
    local out_dir="$RESULT_ROOT/unit_tests"
    mkdir -p "$out_dir"

    record_command unit repo pw_k_serial 0 "$out_dir/pw_k_serial.stdout" "$out_dir/pw_k_serial.stderr" \
        "$BUILD_DIR/source/source_basis/module_pw/test_serial/MODULE_PW_basis_pw_k_serial"
    record_command unit repo pw_mpi_np1 0 "$out_dir/pw_mpi_np1.stdout" "$out_dir/pw_mpi_np1.stderr" \
        "$BUILD_DIR/source/source_basis/module_pw/test/MODULE_PW_pw_test"
    record_command unit repo pw_mpi_np2 0 "$out_dir/pw_mpi_np2.stdout" "$out_dir/pw_mpi_np2.stderr" \
        mpirun --allow-run-as-root -np 2 "$BUILD_DIR/source/source_basis/module_pw/test/MODULE_PW_pw_test"
}

run_benchmarks() {
    local -a cases configs
    IFS=',' read -r -a cases <<< "$CASES_CSV"
    IFS=',' read -r -a configs <<< "$CONFIGS_CSV"

    for case_name in "${cases[@]}"; do
        local case_dir="$SCRIPT_DIR/$case_name"
        local case_result="$RESULT_ROOT/$case_name"
        mkdir -p "$case_result"

        if [[ ! -d "$case_dir" ]]; then
            echo "Missing case directory: $case_dir" | tee -a "$ERROR_FILE"
            [[ "$ALLOW_FAILURES" -eq 0 ]] && exit 1
            continue
        fi

        for config in "${configs[@]}"; do
            local np omp label
            IFS=':' read -r np omp label <<< "$config"

            for rep in $(seq 1 "$REPS"); do
                local run_name="np${np}_omp${omp}_r${rep}"
                local stdout_file="$case_result/${run_name}.stdout"
                local stderr_file="$case_result/${run_name}.stderr"
                local start_marker end_marker start_ns end_ns elapsed
                start_marker=$(now_iso)
                start_ns=$(now_ns)

                echo "===== RUN START $start_marker case=$case_name config=$label np=$np omp=$omp rep=$rep =====" > "$stdout_file"
                echo "===== RUN START $start_marker case=$case_name config=$label np=$np omp=$omp rep=$rep =====" > "$stderr_file"

                (
                    cd "$case_dir" || exit 1
                    export OMP_NUM_THREADS="$omp"
                    if [[ "$np" -eq 1 ]]; then
                        timeout "$TIMEOUT_PER_RUN" "$ABACUS_BIN"
                    else
                        timeout "$TIMEOUT_PER_RUN" mpirun --allow-run-as-root -np "$np" "$ABACUS_BIN"
                    fi
                ) >> "$stdout_file" 2>> "$stderr_file"
                local status=$?

                end_marker=$(now_iso)
                end_ns=$(now_ns)
                echo "===== RUN END $end_marker case=$case_name config=$label np=$np omp=$omp rep=$rep status=$status =====" >> "$stdout_file"
                echo "===== RUN END $end_marker case=$case_name config=$label np=$np omp=$omp rep=$rep status=$status =====" >> "$stderr_file"

                elapsed=$(elapsed_seconds "$start_ns" "$end_ns")
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                    benchmark "$case_name" "$label" "$rep" "np=$np omp=$omp $ABACUS_BIN" \
                    "$start_marker" "$end_marker" "$start_ns" "$end_ns" "$elapsed" "$status" >> "$LOG_FILE"

                echo "[$end_marker] benchmark case=$case_name config=$label rep=$rep status=$status elapsed=${elapsed}s"

                if [[ "$status" -ne 0 ]]; then
                    echo "FAILED benchmark case=$case_name config=$label rep=$rep status=$status" >> "$ERROR_FILE"
                    [[ "$ALLOW_FAILURES" -eq 0 ]] && exit "$status"
                fi

                if ! grep -q "TIME STATISTICS" "$stdout_file"; then
                    echo "NO_TIMER case=$case_name config=$label rep=$rep" >> "$ERROR_FILE"
                fi

                rm -rf "$case_dir/OUT.ABACUS"
            done
        done
    done
}

write_environment() {
    {
        echo "repo_root=$REPO_ROOT"
        echo "build_dir=$BUILD_DIR"
        echo "abacus_bin=$ABACUS_BIN"
        echo "result_root=$RESULT_ROOT"
        echo "cases=$CASES_CSV"
        echo "configs=$CONFIGS_CSV"
        echo "reps=$REPS"
        echo "timeout=$TIMEOUT_PER_RUN"
        echo "date=$(now_iso)"
        git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null | sed 's/^/git_head=/'
        git -C "$REPO_ROOT" status --short 2>/dev/null | sed 's/^/git_status=/' || true
    } > "$RESULT_ROOT/environment.txt"
}

write_environment

if [[ "$DO_BUILD" -eq 1 ]]; then
    run_builds
fi

if [[ "$DO_UNIT" -eq 1 ]]; then
    run_unit_tests
fi

if [[ "$DO_BENCH" -eq 1 ]]; then
    run_benchmarks
fi

if [[ "$DO_PARSE" -eq 1 ]]; then
    if [[ "$DO_BENCH" -eq 1 ]]; then
        python3 "$SCRIPT_DIR/parse_timers.py" "$RESULT_ROOT" > "$SUMMARY_FILE" 2> "$RESULT_ROOT/parse_timers.stderr" || {
            echo "parse_timers.py failed; see $RESULT_ROOT/parse_timers.stderr" | tee -a "$ERROR_FILE"
            [[ "$ALLOW_FAILURES" -eq 0 ]] && exit 1
        }
    else
        echo "Benchmark skipped; no timer summary generated." > "$SUMMARY_FILE"
    fi
fi

echo ""
echo "Benchmark workflow complete."
echo "Results: $RESULT_ROOT"
echo "Log:     $LOG_FILE"
echo "Summary: $SUMMARY_FILE"
if [[ -s "$ERROR_FILE" ]]; then
    echo "Errors:  $ERROR_FILE"
fi
