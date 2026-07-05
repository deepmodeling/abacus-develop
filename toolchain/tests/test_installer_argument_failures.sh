#!/usr/bin/env bash
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
TOOLCHAIN_DIR="${REPO_ROOT}/toolchain"
FAILURES=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  FAILURES=$((FAILURES + 1))
}

copy_toolchain() {
  local tmpdir="$1"
  local entry name
  mkdir -p "${tmpdir}/toolchain"
  while IFS= read -r -d '' entry; do
    name="${entry##*/}"
    case "$name" in
      build|install) continue ;;
    esac
    cp -a "$entry" "${tmpdir}/toolchain/"
  done < <(find "${TOOLCHAIN_DIR}" -mindepth 1 -maxdepth 1 -print0)
}

run_installer_in_copy() {
  local tmpdir="$1"
  shift
  (cd "${tmpdir}/toolchain" && ./install_abacus_toolchain_new.sh "$@") >"${tmpdir}/output.log" 2>&1
}

assert_invalid_input_fails() {
  local name="$1"
  local expected_text="$2"
  shift 2

  local tmpdir status
  tmpdir="$(mktemp -d)"
  copy_toolchain "$tmpdir"

  run_installer_in_copy "$tmpdir" "$@"
  status=$?

  if [[ "$status" -eq 0 ]]; then
    cat "${tmpdir}/output.log" >&2
    fail "${name} exited 0; expected nonzero"
  fi

  if ! grep -Fq -- "$expected_text" "${tmpdir}/output.log"; then
    cat "${tmpdir}/output.log" >&2
    fail "${name} did not report expected error: ${expected_text}"
  fi

  if ! grep -Fq "install_abacus_toolchain_new.sh [OPTIONS]" "${tmpdir}/output.log"; then
    cat "${tmpdir}/output.log" >&2
    fail "${name} output did not contain usage text"
  fi

  if [[ -e "${tmpdir}/toolchain/install/setup" ]]; then
    fail "${name} wrote install/setup even though argument parsing failed"
  fi

  rm -rf "$tmpdir"
}

assert_invalid_input_fails "invalid package version" "Invalid package version format" --dry-run --package-version bad:wrong
assert_invalid_input_fails "invalid gpu version" "Invalid GPU version" --dry-run --gpu-ver bad

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s installer argument failure test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'installer argument failure tests passed\n'
