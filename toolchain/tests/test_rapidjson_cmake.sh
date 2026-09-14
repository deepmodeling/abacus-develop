#!/usr/bin/env bash
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
FAILURES=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  FAILURES=$((FAILURES + 1))
}

write_target_package() {
  local prefix="$1"
  mkdir -p "${prefix}/include/rapidjson" "${prefix}/lib/cmake/RapidJSON"
  printf '#pragma once\n' >"${prefix}/include/rapidjson/document.h"
  cat >"${prefix}/lib/cmake/RapidJSON/RapidJSONConfig.cmake" <<'EOF'
message(STATUS "Loaded fake RapidJSON target package")
get_filename_component(RapidJSON_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_rapidjson_prefix "${RapidJSON_CMAKE_DIR}/../../.." ABSOLUTE)
if(NOT TARGET RapidJSON)
  add_library(RapidJSON INTERFACE IMPORTED)
  set_property(TARGET RapidJSON PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${_rapidjson_prefix}/include")
endif()
EOF
}

write_fake_mkl() {
  local prefix="$1"
  mkdir -p "${prefix}/include" "${prefix}/lib"
  printf '#pragma once\n' >"${prefix}/include/mkl_service.h"
  : >"${prefix}/lib/libmkl_core.so"
  : >"${prefix}/lib/libmkl_gf_lp64.so"
  : >"${prefix}/lib/libmkl_gnu_thread.so"
}

write_variable_only_package() {
  local prefix="$1"
  mkdir -p "${prefix}/include/rapidjson" "${prefix}/lib/cmake/RapidJSON"
  printf '#pragma once\n' >"${prefix}/include/rapidjson/document.h"
  cat >"${prefix}/lib/cmake/RapidJSON/RapidJSONConfig.cmake" <<'EOF'
get_filename_component(RapidJSON_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_rapidjson_prefix "${RapidJSON_CMAKE_DIR}/../../.." ABSOLUTE)
set(RapidJSON_INCLUDE_DIR "${_rapidjson_prefix}/include")
set(RapidJSON_INCLUDE_DIRS "${RapidJSON_INCLUDE_DIR}")
EOF
}

run_top_level_configure() {
  local build_dir="$1"
  local prefix="$2"
  local mkl_root="$3"

  cmake -S "$REPO_ROOT" -B "$build_dir" \
    -DENABLE_RAPIDJSON=ON \
    -DENABLE_LCAO=OFF \
    -DENABLE_MPI=OFF \
    -DENABLE_OPENMP=OFF \
    -DMKLROOT="$mkl_root" \
    -DCMAKE_PREFIX_PATH="$prefix" \
    >"${build_dir}.log" 2>&1
}

test_top_level_accepts_rapidjson_target_package() {
  local tmpdir prefix mkl_root build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  mkl_root="${tmpdir}/mkl"
  build_dir="${tmpdir}/target-build"

  write_target_package "$prefix"
  write_fake_mkl "$mkl_root"
  run_top_level_configure "$build_dir" "$prefix" "$mkl_root"
  status=$?

  if ! grep -Fq "Loaded fake RapidJSON target package" "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake did not load the fake RapidJSON target package"
  elif grep -Fq "RapidJSON was found, but target RapidJSON is missing" "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake rejected RapidJSON target package as target-missing"
  elif grep -Fq 'Could not find a package configuration file provided by "RapidJSON"' "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake did not find the fake RapidJSON target package"
  elif [[ "$status" -ne 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake failed with a RapidJSON target package"
  fi

  rm -rf "$tmpdir"
}

test_top_level_rejects_variable_only_package() {
  local tmpdir prefix mkl_root build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  mkl_root="${tmpdir}/mkl"
  build_dir="${tmpdir}/variable-build"

  write_variable_only_package "$prefix"
  write_fake_mkl "$mkl_root"
  run_top_level_configure "$build_dir" "$prefix" "$mkl_root"
  status=$?

  if [[ "$status" -eq 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake configured with variable-only RapidJSON package; expected target-missing failure"
  elif ! grep -Fq "RapidJSON was found, but target RapidJSON is missing." "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake failed for the wrong reason with variable-only RapidJSON package"
  fi

  rm -rf "$tmpdir"
}

test_top_level_accepts_rapidjson_target_package
test_top_level_rejects_variable_only_package

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s RapidJSON CMake test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'RapidJSON CMake tests passed\n'
