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
get_filename_component(RapidJSON_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_rapidjson_prefix "${RapidJSON_CMAKE_DIR}/../../.." ABSOLUTE)
if(NOT TARGET RapidJSON)
  add_library(RapidJSON INTERFACE IMPORTED)
  set_property(TARGET RapidJSON PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${_rapidjson_prefix}/include")
endif()
EOF
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

write_test_project() {
  local project_dir="$1"
  mkdir -p "$project_dir"
  cat >"${project_dir}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(check_rapidjson_target LANGUAGES CXX)

function(abacus_add_feature_definitions)
  set_property(GLOBAL APPEND PROPERTY ABACUS_FEATURE_DEFINITIONS ${ARGN})
endfunction()

find_package(RapidJSON CONFIG REQUIRED)
if(NOT TARGET RapidJSON)
  message(FATAL_ERROR "RapidJSON package configuration did not define the RapidJSON target")
endif()

add_library(deps INTERFACE)
target_link_libraries(deps INTERFACE RapidJSON)
abacus_add_feature_definitions(__RAPIDJSON)

get_target_property(_links deps INTERFACE_LINK_LIBRARIES)
if(NOT "RapidJSON" IN_LIST _links)
  message(FATAL_ERROR "Expected deps to link to RapidJSON target, got '${_links}'")
endif()

get_property(_defs GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS)
if(NOT "__RAPIDJSON" IN_LIST _defs)
  message(FATAL_ERROR "Expected __RAPIDJSON feature definition")
endif()
EOF
}

run_cmake_configure() {
  local source_dir="$1"
  local build_dir="$2"
  local prefix="$3"

  cmake -S "$source_dir" -B "$build_dir" \
    -DCMAKE_PREFIX_PATH="$prefix" \
    >"${build_dir}.log" 2>&1
}

run_top_level_configure() {
  local build_dir="$1"
  local prefix="$2"

  cmake -S "$REPO_ROOT" -B "$build_dir" \
    -DENABLE_RAPIDJSON=ON \
    -DCMAKE_PREFIX_PATH="$prefix" \
    >"${build_dir}.log" 2>&1
}

test_rapidjson_target_package_configures() {
  local tmpdir prefix source_dir build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  source_dir="${tmpdir}/target-project"
  build_dir="${tmpdir}/target-build"

  write_target_package "$prefix"
  write_test_project "$source_dir"
  run_cmake_configure "$source_dir" "$build_dir" "$prefix"
  status=$?

  if [[ "$status" -ne 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "RapidJSON target package did not configure"
  fi

  rm -rf "$tmpdir"
}

test_variable_only_package_fails() {
  local tmpdir prefix source_dir build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  source_dir="${tmpdir}/variable-project"
  build_dir="${tmpdir}/variable-build"

  write_variable_only_package "$prefix"
  write_test_project "$source_dir"
  run_cmake_configure "$source_dir" "$build_dir" "$prefix"
  status=$?

  if [[ "$status" -eq 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "variable-only RapidJSON package configured; expected target-missing failure"
  elif ! grep -Fq "RapidJSON package configuration did not define the RapidJSON target" "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "variable-only RapidJSON package failed for the wrong reason"
  fi

  rm -rf "$tmpdir"
}

test_top_level_rejects_variable_only_package() {
  local tmpdir prefix build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  build_dir="${tmpdir}/top-level-build"

  write_variable_only_package "$prefix"
  run_top_level_configure "$build_dir" "$prefix"
  status=$?

  if [[ "$status" -eq 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake configured with variable-only RapidJSON package; expected target-missing failure"
  elif ! grep -Fq "RapidJSON package configuration did not define the RapidJSON target" "${build_dir}.log"; then
    cat "${build_dir}.log" >&2
    fail "top-level CMake failed for the wrong reason with variable-only RapidJSON package"
  fi

  rm -rf "$tmpdir"
}

test_top_level_uses_rapidjson_target_directly() {
  local cmakelists="${REPO_ROOT}/CMakeLists.txt"

  if [[ -e "${REPO_ROOT}/cmake/AbacusRapidJSON.cmake" ]]; then
    fail "cmake/AbacusRapidJSON.cmake still exists; RapidJSON config should provide the target"
  fi

  if ! grep -Fq "find_package(RapidJSON CONFIG REQUIRED)" "$cmakelists"; then
    fail "top-level CMakeLists.txt does not call find_package(RapidJSON CONFIG REQUIRED)"
  fi

  if ! grep -Fq "target_link_libraries(abacus_external_deps INTERFACE RapidJSON)" "$cmakelists"; then
    fail "top-level CMakeLists.txt does not link abacus_external_deps to RapidJSON target"
  fi

  if grep -Fq "include(cmake/AbacusRapidJSON.cmake)" "$cmakelists"; then
    fail "top-level CMakeLists.txt still includes the RapidJSON compatibility helper"
  fi
}

test_aocc_build_script_does_not_override_rapidjson_dir() {
  local script="${REPO_ROOT}/toolchain/build_abacus_aocc-aocl.sh"
  if grep -Fq -- '-DRapidJSON_DIR=$RAPIDJSON' "$script"; then
    fail "${script} still passes RapidJSON_DIR as the package prefix"
  fi
}

test_rapidjson_target_package_configures
test_variable_only_package_fails
test_top_level_rejects_variable_only_package
test_top_level_uses_rapidjson_target_directly
test_aocc_build_script_does_not_override_rapidjson_dir

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s RapidJSON CMake test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'RapidJSON CMake tests passed\n'
