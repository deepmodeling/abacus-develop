#!/usr/bin/env bash
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
FAILURES=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  FAILURES=$((FAILURES + 1))
}

write_lowercase_target_package() {
  local prefix="$1"
  mkdir -p "${prefix}/include/rapidjson" "${prefix}/lib/cmake/RapidJSON"
  printf '#pragma once\n' >"${prefix}/include/rapidjson/document.h"
  cat >"${prefix}/lib/cmake/RapidJSON/RapidJSONConfig.cmake" <<'EOF'
get_filename_component(RapidJSON_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_rapidjson_prefix "${RapidJSON_CMAKE_DIR}/../../.." ABSOLUTE)
set(RapidJSON_INCLUDE_DIR "${_rapidjson_prefix}/include")
set(RapidJSON_INCLUDE_DIRS "${RapidJSON_INCLUDE_DIR}")
if(NOT TARGET rapidjson)
  add_library(rapidjson INTERFACE IMPORTED)
  set_property(TARGET rapidjson PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${RapidJSON_INCLUDE_DIRS}")
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

write_target_test_project() {
  local project_dir="$1"
  mkdir -p "$project_dir"
  cat >"${project_dir}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(check_rapidjson_target LANGUAGES CXX)

function(abacus_add_feature_definitions)
  set_property(GLOBAL APPEND PROPERTY ABACUS_FEATURE_DEFINITIONS ${ARGN})
endfunction()

include("${ABACUS_SOURCE_DIR}/cmake/AbacusRapidJSON.cmake")

add_library(deps INTERFACE)
abacus_configure_rapidjson(deps)

get_target_property(_links deps INTERFACE_LINK_LIBRARIES)
if(NOT _links MATCHES "(^|;)rapidjson($|;)")
  message(FATAL_ERROR "Expected deps to link to rapidjson target, got '${_links}'")
endif()

get_property(_defs GLOBAL PROPERTY ABACUS_FEATURE_DEFINITIONS)
if(NOT "__RAPIDJSON" IN_LIST _defs)
  message(FATAL_ERROR "Expected __RAPIDJSON feature definition")
endif()
EOF
}

write_variable_test_project() {
  local project_dir="$1"
  mkdir -p "$project_dir"
  cat >"${project_dir}/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(check_rapidjson_variables LANGUAGES CXX)

function(abacus_add_feature_definitions)
  set_property(GLOBAL APPEND PROPERTY ABACUS_FEATURE_DEFINITIONS ${ARGN})
endfunction()

include("${ABACUS_SOURCE_DIR}/cmake/AbacusRapidJSON.cmake")

add_library(deps INTERFACE)
abacus_configure_rapidjson(deps)

get_target_property(_includes deps INTERFACE_INCLUDE_DIRECTORIES)
if(NOT _includes STREQUAL "${EXPECTED_RAPIDJSON_INCLUDE}")
  message(FATAL_ERROR "Expected include '${EXPECTED_RAPIDJSON_INCLUDE}', got '${_includes}'")
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
  local expected_include="$4"

  cmake -S "$source_dir" -B "$build_dir" \
    -DABACUS_SOURCE_DIR="${REPO_ROOT}" \
    -DCMAKE_PREFIX_PATH="$prefix" \
    -DEXPECTED_RAPIDJSON_INCLUDE="$expected_include" \
    >"${build_dir}.log" 2>&1
}

test_lowercase_target_package() {
  local tmpdir prefix source_dir build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  source_dir="${tmpdir}/target-project"
  build_dir="${tmpdir}/target-build"

  write_lowercase_target_package "$prefix"
  write_target_test_project "$source_dir"
  run_cmake_configure "$source_dir" "$build_dir" "$prefix" "${prefix}/include"
  status=$?

  if [[ "$status" -ne 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "lowercase RapidJSON target package did not configure"
  fi

  rm -rf "$tmpdir"
}

test_variable_only_package() {
  local tmpdir prefix source_dir build_dir status
  tmpdir="$(mktemp -d)"
  prefix="${tmpdir}/prefix"
  source_dir="${tmpdir}/variable-project"
  build_dir="${tmpdir}/variable-build"

  write_variable_only_package "$prefix"
  write_variable_test_project "$source_dir"
  run_cmake_configure "$source_dir" "$build_dir" "$prefix" "${prefix}/include"
  status=$?

  if [[ "$status" -ne 0 ]]; then
    cat "${build_dir}.log" >&2
    fail "variable-only RapidJSON package did not configure"
  fi

  rm -rf "$tmpdir"
}

test_aocc_build_script_does_not_override_rapidjson_dir() {
  local script="${REPO_ROOT}/toolchain/build_abacus_aocc-aocl.sh"
  if grep -Fq -- '-DRapidJSON_DIR=$RAPIDJSON' "$script"; then
    fail "${script} still passes RapidJSON_DIR as the package prefix"
  fi
}

test_lowercase_target_package
test_variable_only_package
test_aocc_build_script_does_not_override_rapidjson_dir

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s RapidJSON CMake test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'RapidJSON CMake tests passed\n'
