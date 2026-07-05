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
  set_property(TARGET RapidJSON PROPERTY INTERFACE_COMPILE_DEFINITIONS ABACUS_TEST_RAPIDJSON_TARGET_USAGE)
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

  mkdir -p "${build_dir}/.cmake/api/v1/query"
  touch "${build_dir}/.cmake/api/v1/query/codemodel-v2"

  cmake -S "$REPO_ROOT" -B "$build_dir" \
    -DENABLE_RAPIDJSON=ON \
    -DENABLE_LCAO=OFF \
    -DENABLE_MPI=OFF \
    -DUSE_OPENMP=OFF \
    -DMKLROOT="$mkl_root" \
    -DCMAKE_PREFIX_PATH="$prefix" \
    >"${build_dir}.log" 2>&1
}

assert_abacus_target_has_rapidjson_usage() {
  local build_dir="$1"

  python3 - "$build_dir" <<'PY'
import json
import pathlib
import sys

build_dir = pathlib.Path(sys.argv[1])
reply_dir = build_dir / ".cmake" / "api" / "v1" / "reply"

try:
    index_file = next(reply_dir.glob("index-*.json"))
    index = json.loads(index_file.read_text())
    codemodel_file = reply_dir / index["reply"]["codemodel-v2"]["jsonFile"]
    codemodel = json.loads(codemodel_file.read_text())
    target_refs = [
        target
        for target in codemodel["configurations"][0]["targets"]
        if target["name"].startswith("abacus_")
    ]
except (KeyError, StopIteration, FileNotFoundError, json.JSONDecodeError) as exc:
    print(f"failed to read CMake File API codemodel: {exc}", file=sys.stderr)
    sys.exit(1)

if len(target_refs) != 1:
    names = ", ".join(sorted(target["name"] for target in target_refs)) or "<none>"
    print(f"expected exactly one ABACUS executable target, got: {names}", file=sys.stderr)
    sys.exit(1)

target_ref = target_refs[0]
target = json.loads((reply_dir / target_ref["jsonFile"]).read_text())
definitions = {
    definition["define"]
    for group in target.get("compileGroups", [])
    for definition in group.get("defines", [])
    if "define" in definition
}
required = {"__RAPIDJSON", "ABACUS_TEST_RAPIDJSON_TARGET_USAGE"}
missing = sorted(required - definitions)

if missing:
    observed = ", ".join(sorted(definitions)) or "<none>"
    print(
        f"{target_ref['name']} is missing RapidJSON definitions: {', '.join(missing)}",
        file=sys.stderr,
    )
    print(f"observed definitions: {observed}", file=sys.stderr)
    sys.exit(1)
PY
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
  elif ! assert_abacus_target_has_rapidjson_usage "$build_dir"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake did not propagate RapidJSON feature and target usage"
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
