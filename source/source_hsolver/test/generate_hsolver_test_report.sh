#!/bin/bash

set -o pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

find_build_dir() {
    local dir="$1"
    while [[ "$dir" != "/" ]]; do
        if [[ -f "$dir/CMakeCache.txt" && -f "$dir/CTestTestfile.cmake" ]]; then
            echo "$dir"
            return 0
        fi
        dir=$(dirname "$dir")
    done
    return 1
}

build_dir=${HSOLVER_TEST_BUILD_DIR:-}
if [[ -z "$build_dir" ]]; then
    build_dir=$(find_build_dir "$PWD")
fi
if [[ -z "$build_dir" ]]; then
    build_dir=$(find_build_dir "$script_dir")
fi
if [[ -z "$build_dir" ]]; then
    echo "Cannot locate the CTest build directory. Set HSOLVER_TEST_BUILD_DIR." >&2
    exit 1
fi

report_dir=${HSOLVER_TEST_REPORT_DIR:-"$build_dir/test_reports/hsolver"}
test_regex=${HSOLVER_TEST_REGEX:-"^MODULE_HSOLVER_(ppcg|bpcg)$"}
timestamp=$(date +"%Y%m%d_%H%M%S")

mkdir -p "$report_dir"

xml_file="$report_dir/hsolver_unit_tests_${timestamp}.xml"
log_file="$report_dir/hsolver_unit_tests_${timestamp}.log"
ctest_has_junit=0
if ctest --help 2>/dev/null | grep -q -- "--output-junit"; then
    ctest_has_junit=1
fi

echo "Build directory : $build_dir"
echo "Test regex      : $test_regex"
if [[ $ctest_has_junit -eq 1 ]]; then
    echo "JUnit XML       : $xml_file"
else
    echo "JUnit XML       : skipped (--output-junit is not supported by this CTest)"
fi
echo "Text log        : $log_file"

if [[ $ctest_has_junit -eq 1 ]]; then
    ctest --test-dir "$build_dir" -V -R "$test_regex" --output-junit "$xml_file" 2>&1 | tee "$log_file"
    status=${PIPESTATUS[0]}
else
    ctest --test-dir "$build_dir" -V -R "$test_regex" 2>&1 | tee "$log_file"
    status=${PIPESTATUS[0]}
fi

if [[ $status -eq 0 ]]; then
    echo "Generated hsolver test reports in: $report_dir"
else
    echo "CTest failed. Partial reports are available in: $report_dir" >&2
fi

exit $status
