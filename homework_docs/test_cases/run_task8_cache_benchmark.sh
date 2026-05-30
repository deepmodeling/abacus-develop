#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_task8_cache_benchmark.sh --case-dir DIR --label LABEL --out-dir DIR [options]

Required:
  --case-dir DIR        ABACUS test case directory.
  --label LABEL         Result label, e.g. baseline or cache.
  --out-dir DIR         Output directory for this benchmark run.

Executable selection:
  --abacus-bin PATH     Path to ABACUS executable.
  --build-dir DIR       Build directory. Used to infer ABACUS executable and compiler info.
  --repo-dir DIR        ABACUS repository root. Defaults to parent of this script's parent.

Run control:
  --nproc N             MPI process count. Default: 1
  --threads N           OpenMP thread count. Default: 1
  --repeat N            Number of measured repeats after warmup. Default: 3
  --warmup N            Number of warmup runs. Default: 1
  --timeout SEC         Timeout per run in seconds. Default: 1800
  --mpi-launch CMD      MPI launcher command. Default: mpirun
  --dry-run             Print resolved configuration and exit.
  -h, --help            Show this help.
EOF
}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
default_repo_dir=$(cd "$script_dir/../.." && pwd)

repo_dir="$default_repo_dir"
build_dir=""
abacus_bin=""
case_dir=""
label=""
out_dir=""
nproc=1
threads=1
repeat=3
warmup=1
timeout_sec=1800
mpi_launch="mpirun"
dry_run=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-dir)
            repo_dir="$2"
            shift 2
            ;;
        --build-dir)
            build_dir="$2"
            shift 2
            ;;
        --abacus-bin)
            abacus_bin="$2"
            shift 2
            ;;
        --case-dir)
            case_dir="$2"
            shift 2
            ;;
        --label)
            label="$2"
            shift 2
            ;;
        --nproc)
            nproc="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --repeat)
            repeat="$2"
            shift 2
            ;;
        --warmup)
            warmup="$2"
            shift 2
            ;;
        --timeout)
            timeout_sec="$2"
            shift 2
            ;;
        --out-dir)
            out_dir="$2"
            shift 2
            ;;
        --mpi-launch)
            mpi_launch="$2"
            shift 2
            ;;
        --dry-run)
            dry_run=1
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

if [[ -z "$case_dir" || -z "$label" || -z "$out_dir" ]]; then
    echo "Missing required arguments." >&2
    usage >&2
    exit 1
fi

if [[ -n "$build_dir" && -n "$abacus_bin" ]]; then
    echo "Please provide either --build-dir or --abacus-bin, not both." >&2
    exit 1
fi

if [[ -z "$abacus_bin" && -z "$build_dir" ]]; then
    echo "Please provide either --abacus-bin or --build-dir." >&2
    exit 1
fi

if [[ -n "$build_dir" ]]; then
    build_dir=$(cd "$build_dir" && pwd)
fi
repo_dir=$(cd "$repo_dir" && pwd)
case_dir=$(cd "$case_dir" && pwd)
mkdir -p "$out_dir"
out_dir=$(cd "$out_dir" && pwd)

resolve_abacus_bin() {
    local from_build_dir="$1"
    local candidate
    for candidate in \
        "$from_build_dir/abacus" \
        "$from_build_dir/abacus_basic_para" \
        "$from_build_dir/ABACUS" \
        "$from_build_dir/bin/abacus" \
        "$from_build_dir/source/abacus"; do
        if [[ -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    candidate=$(find "$from_build_dir" -maxdepth 2 -type f -perm -111 -name 'abacus*' | sort | head -n 1 || true)
    if [[ -n "$candidate" ]]; then
        printf '%s\n' "$candidate"
        return 0
    fi
    return 1
}

if [[ -z "$abacus_bin" ]]; then
    if ! abacus_bin=$(resolve_abacus_bin "$build_dir"); then
        echo "Could not infer ABACUS executable from build dir: $build_dir" >&2
        exit 1
    fi
else
    abacus_bin=$(cd "$(dirname "$abacus_bin")" && pwd)/$(basename "$abacus_bin")
fi

if [[ ! -x "$abacus_bin" ]]; then
    echo "ABACUS executable is not runnable: $abacus_bin" >&2
    exit 1
fi

if [[ ! -d "$case_dir" ]]; then
    echo "Case directory does not exist: $case_dir" >&2
    exit 1
fi

run_csv="$out_dir/run_manifest.csv"
env_txt="$out_dir/environment_info.txt"
run_log="$out_dir/run_log.txt"
logs_dir="$out_dir/logs"
mkdir -p "$logs_dir"

case_name=$(basename "$case_dir")
timestamp=$(date '+%Y-%m-%d %H:%M:%S %z')
hostname_value=$(hostname 2>/dev/null || uname -n)
git_branch=$(git -C "$repo_dir" branch --show-current 2>/dev/null || echo "NA")
git_commit=$(git -C "$repo_dir" rev-parse HEAD 2>/dev/null || echo "NA")
git_status_short=$(git -C "$repo_dir" status --short 2>/dev/null || true)

compiler_c="NA"
compiler_cxx="NA"
compiler_c_id="NA"
compiler_cxx_id="NA"
compiler_cxx_ver="NA"
if [[ -n "$build_dir" && -f "$build_dir/CMakeCache.txt" ]]; then
    compiler_c=$(grep -E '^CMAKE_C_COMPILER:FILEPATH=' "$build_dir/CMakeCache.txt" | head -n1 | cut -d= -f2- || true)
    compiler_cxx=$(grep -E '^CMAKE_CXX_COMPILER:FILEPATH=' "$build_dir/CMakeCache.txt" | head -n1 | cut -d= -f2- || true)
    compiler_c_id=$(grep -E '^CMAKE_C_COMPILER_ID:STRING=' "$build_dir/CMakeCache.txt" | head -n1 | cut -d= -f2- || true)
    compiler_cxx_id=$(grep -E '^CMAKE_CXX_COMPILER_ID:STRING=' "$build_dir/CMakeCache.txt" | head -n1 | cut -d= -f2- || true)
    compiler_cxx_ver=$(grep -E '^CMAKE_CXX_COMPILER_VERSION:STRING=' "$build_dir/CMakeCache.txt" | head -n1 | cut -d= -f2- || true)
fi
mpi_version=$("$mpi_launch" --version 2>&1 | head -n 5 || true)

{
    echo "timestamp=$timestamp"
    echo "hostname=$hostname_value"
    echo "repo_dir=$repo_dir"
    echo "build_dir=${build_dir:-NA}"
    echo "abacus_bin=$abacus_bin"
    echo "case_dir=$case_dir"
    echo "case_name=$case_name"
    echo "label=$label"
    echo "nproc=$nproc"
    echo "threads=$threads"
    echo "warmup=$warmup"
    echo "repeat=$repeat"
    echo "timeout_sec=$timeout_sec"
    echo "mpi_launch=$mpi_launch"
    echo "git_branch=$git_branch"
    echo "git_commit=$git_commit"
    echo "omp_num_threads=$threads"
    echo "omp_proc_bind=${OMP_PROC_BIND:-NA}"
    echo "omp_places=${OMP_PLACES:-NA}"
    echo "omp_schedule=${OMP_SCHEDULE:-NA}"
    echo "omp_dynamic=${OMP_DYNAMIC:-NA}"
    echo "omp_display_env=${OMP_DISPLAY_ENV:-NA}"
    echo "compiler_c=${compiler_c:-NA}"
    echo "compiler_c_id=${compiler_c_id:-NA}"
    echo "compiler_cxx=${compiler_cxx:-NA}"
    echo "compiler_cxx_id=${compiler_cxx_id:-NA}"
    echo "compiler_cxx_version=${compiler_cxx_ver:-NA}"
    echo "uname=$(uname -a)"
    echo "lscpu=$(lscpu 2>/dev/null | tr '\n' ';' || echo NA)"
    echo "mpirun_version=$(printf '%s' "$mpi_version" | tr '\n' ';')"
    echo "git_status_short<<EOF"
    printf '%s\n' "${git_status_short:-}"
    echo "EOF"
} > "$env_txt"

printf 'label,case_name,git_branch,git_commit,nproc,threads,run_kind,run_index,stdout_path,stderr_path,out_abacus_path,exit_code,start_time,end_time\n' > "$run_csv"

{
    echo "Task8 cache benchmark start: $timestamp"
    echo "label=$label case=$case_name nproc=$nproc threads=$threads warmup=$warmup repeat=$repeat"
    echo "abacus_bin=$abacus_bin"
    echo "case_dir=$case_dir"
    echo "out_dir=$out_dir"
} > "$run_log"

if [[ "$dry_run" -eq 1 ]]; then
    cat "$env_txt"
    exit 0
fi

run_one() {
    local run_kind="$1"
    local run_index="$2"
    local run_id="${label}_${case_name}_np${nproc}_omp${threads}_${run_kind}${run_index}"
    local stdout_file="$logs_dir/${run_id}.stdout"
    local stderr_file="$logs_dir/${run_id}.stderr"
    local out_abacus_dir="$logs_dir/${run_id}.OUT.ABACUS"
    local start_time end_time exit_code

    start_time=$(date '+%Y-%m-%d %H:%M:%S %z')
    echo "[$start_time] starting $run_id" | tee -a "$run_log"

    rm -rf "$out_abacus_dir"

    pushd "$case_dir" >/dev/null
    rm -rf OUT.ABACUS
    export OMP_NUM_THREADS="$threads"
    set +e
    if [[ "$nproc" -eq 1 ]]; then
        timeout "$timeout_sec" "$abacus_bin" >"$stdout_file" 2>"$stderr_file"
        exit_code=$?
    else
        timeout "$timeout_sec" "$mpi_launch" -np "$nproc" "$abacus_bin" >"$stdout_file" 2>"$stderr_file"
        exit_code=$?
    fi
    if [[ -d OUT.ABACUS ]]; then
        cp -a OUT.ABACUS "$out_abacus_dir"
    elif [[ -f OUT.ABACUS ]]; then
        cp -a OUT.ABACUS "$out_abacus_dir"
    else
        out_abacus_dir="NA"
    fi
    rm -rf OUT.ABACUS
    set -e
    popd >/dev/null

    end_time=$(date '+%Y-%m-%d %H:%M:%S %z')
    echo "[$end_time] finished $run_id exit_code=$exit_code" | tee -a "$run_log"

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$label" "$case_name" "$git_branch" "$git_commit" "$nproc" "$threads" \
        "$run_kind" "$run_index" "$stdout_file" "$stderr_file" "$out_abacus_dir" "$exit_code" \
        "$start_time" "$end_time" >> "$run_csv"
}

for i in $(seq 1 "$warmup"); do
    run_one "warmup" "$i"
done

for i in $(seq 1 "$repeat"); do
    run_one "repeat" "$i"
done

python3 "$script_dir/collect_task8_cache_results.py" \
    --input-dir "$out_dir" \
    --out-dir "$out_dir"

echo "Benchmark finished. Outputs are under: $out_dir"
