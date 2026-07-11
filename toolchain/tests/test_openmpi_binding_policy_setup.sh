#!/usr/bin/env bash
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
TOOLCHAIN_DIR="${REPO_ROOT}/toolchain"
FAILURES=0

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  FAILURES=$((FAILURES + 1))
}

write_fake_openmpi_commands() {
  local bindir="$1"
  local version="$2"

  mkdir -p "$bindir"
  cat >"${bindir}/mpiexec" <<EOF
#!/usr/bin/env bash
if [[ "\$1" == "--version" ]]; then
  printf 'mpiexec (Open MPI) ${version}\n'
  exit 0
fi
exit 0
EOF
  cat >"${bindir}/mpicxx" <<'EOF'
#!/usr/bin/env bash
if [[ "$1" == "--showme:libs" ]]; then
  printf 'mpi\n'
  exit 0
fi
exit 0
EOF
  cp "${bindir}/mpicxx" "${bindir}/mpicc"
  cp "${bindir}/mpicxx" "${bindir}/mpifort"
  chmod +x "${bindir}/mpiexec" "${bindir}/mpicc" "${bindir}/mpicxx" "${bindir}/mpifort"
}

run_openmpi_system_setup() {
  local tmpdir="$1"
  local version="$2"

  mkdir -p "${tmpdir}/install" "${tmpdir}/build"
  : >"${tmpdir}/install/setup"
  : >"${tmpdir}/install/toolchain.env"
  cat >"${tmpdir}/install/toolchain.conf" <<'EOF'
MPI_MODE="openmpi"
with_openmpi="__SYSTEM__"
PACK_RUN="__FALSE__"
EOF

  write_fake_openmpi_commands "${tmpdir}/fake-bin" "$version"
  PATH="${tmpdir}/fake-bin:${PATH}" \
    ROOTDIR="$tmpdir" \
    SCRIPTDIR="${TOOLCHAIN_DIR}/scripts" \
    INSTALLDIR="${tmpdir}/install" \
    BUILDDIR="${tmpdir}/build" \
    SETUPFILE="${tmpdir}/install/setup" \
    bash "${TOOLCHAIN_DIR}/scripts/stage1/install_openmpi.sh" \
    >"${tmpdir}/openmpi.log" 2>&1
}

assert_setup_contains() {
  local file="$1"
  local expected="$2"

  if ! grep -Fq "$expected" "$file"; then
    cat "$file" >&2
    fail "${file} does not contain expected text: ${expected}"
  fi
}

assert_setup_not_contains() {
  local file="$1"
  local unexpected="$2"

  if grep -Fq "$unexpected" "$file"; then
    cat "$file" >&2
    fail "${file} contains unexpected text: ${unexpected}"
  fi
}

test_openmpi5_setup_disables_prte_binding() {
  local tmpdir status
  tmpdir="$(mktemp -d)"

  run_openmpi_system_setup "$tmpdir" "5.0.10"
  status=$?

  if [[ "$status" -ne 0 ]]; then
    cat "${tmpdir}/openmpi.log" >&2
    fail "OpenMPI 5 setup generation failed with status ${status}"
  else
    assert_setup_contains "${tmpdir}/install/setup" "export PRTE_MCA_hwloc_default_binding_policy=none"
    assert_setup_not_contains "${tmpdir}/install/setup" "export OMPI_MCA_hwloc_base_binding_policy=none"
  fi

  rm -rf "$tmpdir"
}

test_openmpi4_setup_disables_ompi_binding() {
  local tmpdir status
  tmpdir="$(mktemp -d)"

  run_openmpi_system_setup "$tmpdir" "4.1.8"
  status=$?

  if [[ "$status" -ne 0 ]]; then
    cat "${tmpdir}/openmpi.log" >&2
    fail "OpenMPI 4 setup generation failed with status ${status}"
  else
    assert_setup_contains "${tmpdir}/install/setup" "export OMPI_MCA_hwloc_base_binding_policy=none"
    assert_setup_not_contains "${tmpdir}/install/setup" "export PRTE_MCA_hwloc_default_binding_policy=none"
  fi

  rm -rf "$tmpdir"
}

test_openmpi5_setup_disables_prte_binding
test_openmpi4_setup_disables_ompi_binding

if [[ "$FAILURES" -ne 0 ]]; then
  printf '%s OpenMPI binding policy setup test(s) failed\n' "$FAILURES" >&2
  exit 1
fi

printf 'OpenMPI binding policy setup tests passed\n'
