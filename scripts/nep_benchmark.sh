#!/usr/bin/env bash
# NEP parallel benchmark and correctness check for task 6.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NEP_DIR="${NEP_DIR:-${ROOT}/third_party/NEP_CPU/build}"
NEP_MODEL="${NEP_MODEL:-${ROOT}/tests/PP_ORB/nep_hfo2.txt}"
BENCHMARK_BIN="${BENCHMARK_BIN:-${ROOT}/tools/nep_benchmark/build/nep_benchmark}"

build_nep() {
  echo "==> Building OpenMP-enabled NEP from third_party/NEP_CPU"
  cmake -S "${ROOT}/third_party/NEP_CPU" -B "${ROOT}/third_party/NEP_CPU/build" \
    -DUSE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release
  cmake --build "${ROOT}/third_party/NEP_CPU/build" -j"$(nproc)"
}

build_benchmark() {
  echo "==> Building nep_benchmark tool"
  cmake -S "${ROOT}/tools/nep_benchmark" -B "${ROOT}/tools/nep_benchmark/build" \
    -DNEP_DIR="${ROOT}/third_party/NEP_CPU/build" \
    -DUSE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release
  cmake --build "${ROOT}/tools/nep_benchmark/build" -j"$(nproc)"
}

if [[ ! -f "${NEP_DIR}/libnep.so" && ! -f "${NEP_DIR}/lib/libnep.so" ]]; then
  build_nep
fi

if [[ ! -x "${BENCHMARK_BIN}" ]]; then
  build_benchmark
  BENCHMARK_BIN="${ROOT}/tools/nep_benchmark/build/nep_benchmark"
fi

export LD_LIBRARY_PATH="${ROOT}/third_party/NEP_CPU/build:${LD_LIBRARY_PATH:-}"

echo "==> Correctness: compare 1-thread vs N-thread results (multi-element HfO2 model)"
for nthread in 1 2 4 8; do
  export OMP_NUM_THREADS="${nthread}"
  echo "--- OMP_NUM_THREADS=${nthread}"
  "${BENCHMARK_BIN}" --model "${NEP_MODEL}" --natom 64 --repeat 5 --verify
done

echo "==> Performance sweep"
for nthread in 1 2 4 8; do
  export OMP_NUM_THREADS="${nthread}"
  "${BENCHMARK_BIN}" --model "${NEP_MODEL}" --natom 256 --repeat 20 --perf
done

echo "Done."
