#!/usr/bin/env bash
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
TOOLCHAIN_DIR="${REPO_ROOT}/toolchain"
FAILURES=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  FAILURES=$((FAILURES + 1))
}

assert_file_contains() {
  local file="$1"
  local text="$2"
  if ! grep -Fq "$text" "$file"; then
    fail "${file} does not contain expected text: ${text}"
  fi
}

test_runner_preserves_command_failure() {
  # shellcheck source=/dev/null
  source "${TOOLCHAIN_DIR}/scripts/lib/wrapper_runner.sh"

  local tmpdir log status
  tmpdir="$(mktemp -d)"
  log="${tmpdir}/compile.log"

  run_toolchain_with_log "$log" bash -c 'printf "installer stdout\n"; exit 37'
  status=$?

  if [[ "$status" -ne 37 ]]; then
    fail "run_toolchain_with_log returned ${status}; expected 37"
  fi
  assert_file_contains "$log" "installer stdout"

  rm -rf "$tmpdir"
}

test_wrappers_use_runner() {
  local wrappers=(
    "${TOOLCHAIN_DIR}/toolchain_gnu.sh"
    "${TOOLCHAIN_DIR}/toolchain_intel.sh"
    "${TOOLCHAIN_DIR}/toolchain_gcc-mkl.sh"
    "${TOOLCHAIN_DIR}/toolchain_gcc-aocl.sh"
    "${TOOLCHAIN_DIR}/toolchain_aocc-aocl.sh"
  )

  local wrapper
  for wrapper in "${wrappers[@]}"; do
    assert_file_contains "$wrapper" 'source "${SCRIPT_DIR}/scripts/lib/wrapper_runner.sh"'
    assert_file_contains "$wrapper" 'run_toolchain_with_log compile.log ./install_abacus_toolchain_new.sh'

    if grep -Eq '\|\s*tee[[:space:]]+compile\.log' "$wrapper"; then
      fail "${wrapper} still contains a raw pipe to tee compile.log"
    fi

    if grep -Eq '^exec[[:space:]]+\./install_abacus_toolchain_new\.sh' "$wrapper"; then
      fail "${wrapper} still execs the installer directly"
    fi
  done
}

test_runner_preserves_command_failure
test_wrappers_use_runner

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s wrapper failure propagation test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'wrapper failure propagation tests passed\n'
