#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEFAULT_REPO_DIR=$(cd "$SCRIPT_DIR/../.." && pwd)

# Fixed Task 8 suite configuration.
# Default cases favor controllable runtime while still exercising PW setup,
# plane-wave related data construction, FFT transforms, and SCF iterations.
CASES=(
    "homework_docs/test_cases/gaas_small"
    "homework_docs/test_cases/gaas_medium"
    # "homework_docs/test_cases/gaas_large"
)

NPROCS=(1 2 4)
THREADS=(1 2 4)
WARMUP=1
REPEAT=3
TIMEOUT=1800
CONTINUE_ON_ERROR=1

usage() {
    cat <<'EOF'
Usage:
  run_task8_cache_suite.sh --label LABEL [--build-dir DIR | --abacus-bin PATH] [options]

Required:
  --label LABEL         Result label, usually baseline or cache.

Executable selection:
  --build-dir DIR       Build directory passed through to run_task8_cache_benchmark.sh.
  --abacus-bin PATH     ABACUS executable passed through to run_task8_cache_benchmark.sh.

Optional:
  --repo-dir DIR        Repository root. Default: auto-detect from script location.
  --out-dir DIR         Suite output directory.
                        Default: homework_docs/test_cases/task8_cache_suite_runs/$label_$timestamp
  --mpi-launch CMD      MPI launcher. Default: mpirun
  --dry-run             Print resolved suite matrix without executing.
  -h, --help            Show this help.
EOF
}

REPO_DIR="$DEFAULT_REPO_DIR"
BUILD_DIR=""
ABACUS_BIN=""
LABEL=""
OUT_DIR=""
MPI_LAUNCH="mpirun"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-dir)
            REPO_DIR="$2"
            shift 2
            ;;
        --build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        --abacus-bin)
            ABACUS_BIN="$2"
            shift 2
            ;;
        --label)
            LABEL="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --mpi-launch)
            MPI_LAUNCH="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$LABEL" ]]; then
    echo "Missing required argument: --label" >&2
    usage >&2
    exit 1
fi

if [[ -n "$BUILD_DIR" && -n "$ABACUS_BIN" ]]; then
    echo "Please provide either --build-dir or --abacus-bin, not both." >&2
    exit 1
fi

if [[ -z "$BUILD_DIR" && -z "$ABACUS_BIN" ]]; then
    echo "Please provide either --build-dir or --abacus-bin." >&2
    exit 1
fi

REPO_DIR=$(cd "$REPO_DIR" && pwd)
if [[ -n "$BUILD_DIR" ]]; then
    BUILD_DIR=$(cd "$BUILD_DIR" && pwd)
fi
if [[ -n "$ABACUS_BIN" ]]; then
    ABACUS_BIN=$(cd "$(dirname "$ABACUS_BIN")" && pwd)/$(basename "$ABACUS_BIN")
fi

timestamp=$(date '+%Y%m%d_%H%M%S')
if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="$REPO_DIR/homework_docs/test_cases/task8_cache_suite_runs/${LABEL}_${timestamp}"
fi
mkdir -p "$OUT_DIR"
OUT_DIR=$(cd "$OUT_DIR" && pwd)

SUITE_MANIFEST="$OUT_DIR/suite_manifest.csv"
SUITE_RUN_LOG="$OUT_DIR/suite_run_log.txt"

printf 'label,case_name,case_dir,nproc,threads,out_dir,exit_code,start_time,end_time\n' > "$SUITE_MANIFEST"

{
    echo "Task8 cache suite start: $(date '+%Y-%m-%d %H:%M:%S %z')"
    echo "repo_dir=$REPO_DIR"
    echo "build_dir=${BUILD_DIR:-NA}"
    echo "abacus_bin=${ABACUS_BIN:-NA}"
    echo "label=$LABEL"
    echo "mpi_launch=$MPI_LAUNCH"
    echo "warmup=$WARMUP repeat=$REPEAT timeout=$TIMEOUT continue_on_error=$CONTINUE_ON_ERROR"
    echo "cases=${CASES[*]}"
    echo "nprocs=${NPROCS[*]}"
    echo "threads=${THREADS[*]}"
    echo
} > "$SUITE_RUN_LOG"

print_matrix() {
    echo "Suite configuration:"
    echo "  repo_dir=$REPO_DIR"
    echo "  build_dir=${BUILD_DIR:-NA}"
    echo "  abacus_bin=${ABACUS_BIN:-NA}"
    echo "  label=$LABEL"
    echo "  out_dir=$OUT_DIR"
    echo "  mpi_launch=$MPI_LAUNCH"
    echo "  warmup=$WARMUP repeat=$REPEAT timeout=$TIMEOUT continue_on_error=$CONTINUE_ON_ERROR"
    echo "  cases:"
    for case_rel in "${CASES[@]}"; do
        echo "    - $case_rel"
    done
    echo "  nprocs: ${NPROCS[*]}"
    echo "  threads: ${THREADS[*]}"
    echo "Planned sub-tests:"
    for case_rel in "${CASES[@]}"; do
        local_case_dir="$REPO_DIR/$case_rel"
        case_name=$(basename "$local_case_dir")
        for nproc in "${NPROCS[@]}"; do
            for thread in "${THREADS[@]}"; do
                echo "  - case=$case_name nproc=$nproc threads=$thread out_dir=$OUT_DIR/${case_name}_np${nproc}_omp${thread}"
            done
        done
    done
}

for case_rel in "${CASES[@]}"; do
    if [[ ! -d "$REPO_DIR/$case_rel" ]]; then
        echo "Configured case directory does not exist: $REPO_DIR/$case_rel" >&2
        exit 1
    fi
done

if [[ "$DRY_RUN" -eq 1 ]]; then
    print_matrix
    exit 0
fi

for case_rel in "${CASES[@]}"; do
    case_dir="$REPO_DIR/$case_rel"
    case_name=$(basename "$case_dir")

    for nproc in "${NPROCS[@]}"; do
        for thread in "${THREADS[@]}"; do
            sub_out_dir="$OUT_DIR/${case_name}_np${nproc}_omp${thread}"
            start_time=$(date '+%Y-%m-%d %H:%M:%S %z')

            echo "[$start_time] starting case=$case_name nproc=$nproc threads=$thread" | tee -a "$SUITE_RUN_LOG"

            cmd=(
                bash "$SCRIPT_DIR/run_task8_cache_benchmark.sh"
                --repo-dir "$REPO_DIR"
                --case-dir "$case_dir"
                --label "$LABEL"
                --nproc "$nproc"
                --threads "$thread"
                --warmup "$WARMUP"
                --repeat "$REPEAT"
                --timeout "$TIMEOUT"
                --mpi-launch "$MPI_LAUNCH"
                --out-dir "$sub_out_dir"
            )
            if [[ -n "$BUILD_DIR" ]]; then
                cmd+=(--build-dir "$BUILD_DIR")
            else
                cmd+=(--abacus-bin "$ABACUS_BIN")
            fi

            set +e
            "${cmd[@]}"
            exit_code=$?
            set -e

            end_time=$(date '+%Y-%m-%d %H:%M:%S %z')
            printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
                "$LABEL" "$case_name" "$case_dir" "$nproc" "$thread" "$sub_out_dir" \
                "$exit_code" "$start_time" "$end_time" >> "$SUITE_MANIFEST"

            echo "[$end_time] finished case=$case_name nproc=$nproc threads=$thread exit_code=$exit_code out_dir=$sub_out_dir" | tee -a "$SUITE_RUN_LOG"

            if [[ "$exit_code" -ne 0 && "$CONTINUE_ON_ERROR" -ne 1 ]]; then
                echo "Stopping because a sub-test failed and CONTINUE_ON_ERROR=$CONTINUE_ON_ERROR" | tee -a "$SUITE_RUN_LOG"
                exit "$exit_code"
            fi
        done
    done
done

python3 "$SCRIPT_DIR/collect_task8_cache_suite_results.py" \
    --suite-dir "$OUT_DIR" \
    --out-dir "$OUT_DIR"

echo "Task8 cache suite finished. Results are under: $OUT_DIR"
