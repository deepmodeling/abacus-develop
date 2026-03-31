#!/usr/bin/env bash
set -euo pipefail

# Build ABACUS in two modes:
# 1) no-tests mode (default): preserve original build flow without GTest
# 2) with-tests mode: build and use local/offline GoogleTest
# Usage:
#   bash tools/build_abacus_with_local_gtest.sh --no-tests [abacus_root]
#   bash tools/build_abacus_with_local_gtest.sh --with-tests [abacus_root] /path/to/googletest
#   GTEST_SRC=/path/to/googletest bash tools/build_abacus_with_local_gtest.sh --with-tests [abacus_root]
#   CC=mpicc CXX=mpic++ bash tools/build_abacus_with_local_gtest.sh --with-tests /path/to/googletest
#   bash tools/build_abacus_with_local_gtest.sh --with-tests --cc mpicc --cxx mpic++ /path/to/googletest
#   bash tools/build_abacus_with_local_gtest.sh --with-tests --install /path/to/googletest
#   bash tools/build_abacus_with_local_gtest.sh --with-tests --use-elpa /path/to/googletest
# Defaults are aligned with your known-good build command and can be overridden
# by environment variables: USE_OPENMP_FLAG, ENABLE_LCAO_FLAG, USE_SW_FLAG,
# SW_MATH_PATH, SW_FFT_PATH, CEREAL_INCLUDE_DIR, CXX_FLAGS_EXTRA.

# Recommended usage for SW platform:
# CC=mpicc CXX=mpic++ bash install.sh --with-tests 
# /path/to/googletest/googletest/

ABACUS_ROOT="$PWD"
GTEST_SRC="${GTEST_SRC:-}"
ENABLE_TESTING="OFF"
C_COMPILER="${CC:-}"
CXX_COMPILER="${CXX:-}"
DO_INSTALL="OFF"
USE_ELPA_FLAG="OFF"
USE_OPENMP_FLAG="${USE_OPENMP_FLAG:-OFF}"
ENABLE_LCAO_FLAG="${ENABLE_LCAO_FLAG:-ON}"
USE_SW_FLAG="${USE_SW_FLAG:-ON}"
SW_MATH_PATH="${SW_MATH_PATH:-/usr/sw/yyzlib/xMath-SACA}"
SW_FFT_PATH="${SW_FFT_PATH:-/usr/sw/yyzlib/fftw-3.3.8}"
CEREAL_INCLUDE_DIR="${CEREAL_INCLUDE_DIR:-/home/export/online1/mdt00/shisuan/swhnu/liu/abacus1/3.9.0.19/abacus-develop-3.9.0.19/cereal/include}"
CXX_FLAGS_EXTRA="${CXX_FLAGS_EXTRA:--I/usr/sw/yyzlib/xMath-SACA/include}"

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-tests)
      ENABLE_TESTING="ON"
      shift
      ;;
    --no-tests)
      ENABLE_TESTING="OFF"
      shift
      ;;
    --cc)
      C_COMPILER="${2:-}"
      shift 2
      ;;
    --cxx)
      CXX_COMPILER="${2:-}"
      shift 2
      ;;
    --install)
      DO_INSTALL="ON"
      shift
      ;;
    --no-install)
      DO_INSTALL="OFF"
      shift
      ;;
    --use-elpa)
      USE_ELPA_FLAG="ON"
      shift
      ;;
    --no-elpa)
      USE_ELPA_FLAG="OFF"
      shift
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ${#POSITIONAL_ARGS[@]} -ge 2 ]]; then
  ABACUS_ROOT="${POSITIONAL_ARGS[0]}"
  GTEST_SRC="${POSITIONAL_ARGS[1]}"
elif [[ ${#POSITIONAL_ARGS[@]} -eq 1 ]]; then
  # Single-argument mode:
  # - with-tests: argument is googletest path
  # - no-tests: argument is ABACUS root
  if [[ "${ENABLE_TESTING}" == "ON" ]]; then
    GTEST_SRC="${POSITIONAL_ARGS[0]}"
  else
    ABACUS_ROOT="${POSITIONAL_ARGS[0]}"
  fi
fi

CMAKE_COMPILER_ARGS=()
if [[ -n "${C_COMPILER}" ]]; then
  CMAKE_COMPILER_ARGS+=("-DCMAKE_C_COMPILER=${C_COMPILER}")
fi
if [[ -n "${CXX_COMPILER}" ]]; then
  CMAKE_COMPILER_ARGS+=("-DCMAKE_CXX_COMPILER=${CXX_COMPILER}")
fi

ABACUS_CMAKE_ARGS=()
ABACUS_CMAKE_ARGS+=("-DUSE_ELPA=${USE_ELPA_FLAG}")
ABACUS_CMAKE_ARGS+=("-DUSE_OPENMP=${USE_OPENMP_FLAG}")
ABACUS_CMAKE_ARGS+=("-DENABLE_LCAO=${ENABLE_LCAO_FLAG}")
ABACUS_CMAKE_ARGS+=("-DUSE_SW=${USE_SW_FLAG}")
ABACUS_CMAKE_ARGS+=("-DSW_MATH=${SW_MATH_PATH}")
ABACUS_CMAKE_ARGS+=("-DSW_FFT=${SW_FFT_PATH}")
ABACUS_CMAKE_ARGS+=("-DCEREAL_INCLUDE_DIR=${CEREAL_INCLUDE_DIR}")

ABACUS_CMAKE_ARGS+=("-DCMAKE_CXX_FLAGS=${CXX_FLAGS_EXTRA}")

if [[ ! -f "${ABACUS_ROOT}/CMakeLists.txt" ]]; then
  echo "Error: invalid ABACUS root: ${ABACUS_ROOT}"
  exit 1
fi

NPROC="$(command -v nproc >/dev/null 2>&1 && nproc || getconf _NPROCESSORS_ONLN || echo 4)"

GTEST_PREFIX="${ABACUS_ROOT}/_local/gtest"
GTEST_BUILD_DIR="${ABACUS_ROOT}/_build_gtest"
ABACUS_BUILD_DIR="${ABACUS_ROOT}/build_local_gtest"
ABACUS_INSTALL_PREFIX="${ABACUS_BUILD_DIR}/install"

if [[ "${ENABLE_TESTING}" == "ON" ]]; then
  if [[ -z "${GTEST_SRC}" ]]; then
    echo "Error: with-tests mode requires googletest source path (arg or GTEST_SRC)"
    exit 1
  fi
  if [[ ! -f "${GTEST_SRC}/CMakeLists.txt" ]]; then
    echo "Error: invalid googletest source path: ${GTEST_SRC}"
    exit 1
  fi

  # If user passed repo subdir .../googletest, normalize to repo root so gmock is built too.
  if [[ -d "${GTEST_SRC}/../googlemock" && -f "${GTEST_SRC}/../CMakeLists.txt" ]]; then
    GTEST_SRC="$(cd "${GTEST_SRC}/.." && pwd)"
  fi

  echo "== [1/5] Build and install googletest to ${GTEST_PREFIX}"
  cmake -S "${GTEST_SRC}" -B "${GTEST_BUILD_DIR}" \
    "${CMAKE_COMPILER_ARGS[@]}" \
    -DBUILD_GTEST=ON \
    -DBUILD_GMOCK=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${GTEST_PREFIX}"

  cmake --build "${GTEST_BUILD_DIR}" -j"${NPROC}"

  # Always install googletest to a user-writable local prefix.
  cmake --install "${GTEST_BUILD_DIR}"

  # Use install-tree package config for stable imported targets.
  GTEST_CMAKE_DIR=""
  for d in "${GTEST_PREFIX}/lib/cmake/GTest" "${GTEST_PREFIX}/lib64/cmake/GTest"; do
    if [[ -f "${d}/GTestConfig.cmake" ]]; then
      GTEST_CMAKE_DIR="${d}"
      break
    fi
  done

  if [[ -z "${GTEST_CMAKE_DIR}" ]]; then
    echo "Error: GTestConfig.cmake not found in build or install tree"
    exit 1
  fi

  echo "== [2/5] Configure ABACUS with local GTest (offline mode)"
  cmake -S "${ABACUS_ROOT}" -B "${ABACUS_BUILD_DIR}" \
    "${CMAKE_COMPILER_ARGS[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${ABACUS_INSTALL_PREFIX}" \
    "${ABACUS_CMAKE_ARGS[@]}" \
    -DBUILD_TESTING=ON \
    -DGTest_DIR="${GTEST_CMAKE_DIR}" \
    -DGTEST_DIR="${GTEST_CMAKE_DIR}" \
    -DCMAKE_PREFIX_PATH="${GTEST_PREFIX}" \
    -DFETCHCONTENT_FULLY_DISCONNECTED=ON
else
  echo "== [1/3] Configure ABACUS without tests (original mode)"
  cmake -S "${ABACUS_ROOT}" -B "${ABACUS_BUILD_DIR}" \
    "${CMAKE_COMPILER_ARGS[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${ABACUS_INSTALL_PREFIX}" \
    "${ABACUS_CMAKE_ARGS[@]}" \
    -DBUILD_TESTING=OFF
fi

echo "== Build ABACUS"
cmake --build "${ABACUS_BUILD_DIR}" -j"${NPROC}"

if [[ "${DO_INSTALL}" == "ON" ]]; then
  echo "== Install ABACUS"
  cmake --install "${ABACUS_BUILD_DIR}"
else
  echo "== Skip install (use --install to enable)"
fi

if [[ "${ENABLE_TESTING}" == "ON" ]]; then
  echo "== [5/5] List tests"
  cd "${ABACUS_BUILD_DIR}"
  ctest -N || true
fi

if [[ "${ENABLE_TESTING}" == "ON" ]]; then
  echo "Done. ABACUS is configured to use local GTest at: ${GTEST_CMAKE_DIR}"
else
  echo "Done. ABACUS is built in original mode without GTest (BUILD_TESTING=OFF)."
fi

