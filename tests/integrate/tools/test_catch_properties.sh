#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
integrate_dir=$(realpath "${script_dir}/..")
test_root=$(mktemp -d "${TMPDIR:-/tmp}/catch-properties-test.XXXXXX")
trap 'rm -rf "${test_root}"' EXIT

case_dir="${test_root}/tests/01_PW/test_case"
mkdir -p "${case_dir}/OUT.autotest"
ln -s "${integrate_dir}" "${test_root}/tests/integrate"

cat > "${case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
out_band 1
EOF

cat > "${case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${case_dir}/sr_nao.csr.ref"
touch "${case_dir}/OUT.autotest/sr_nao.csr"

(
    cd "${case_dir}"
    bash -e ../../integrate/tools/catch_properties.sh result.out
)

grep -qx 'CompareBand_pass 1' "${case_dir}/result.out"

if (
    cd "${case_dir}"
    bash ../../integrate/tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a failed comparison while generating result.ref" >&2
    exit 1
fi

if [ -e "${case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial result.ref after comparison failure" >&2
    exit 1
fi

npz_case_dir="${test_root}/tests/01_PW/npz_case"
mkdir -p "${npz_case_dir}/OUT.autotest"

cat > "${npz_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
out_hsr_npz 1
EOF

cat > "${npz_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${npz_case_dir}/sr_nao.csr.ref"
touch "${npz_case_dir}/OUT.autotest/sr_nao.csr"

(
    cd "${npz_case_dir}"
    bash -e ../../integrate/tools/catch_properties.sh result.out
)

grep -qx 'OutputSRNPZ_pass 1' "${npz_case_dir}/result.out"
grep -qx 'OutputHRNPZ_pass 1' "${npz_case_dir}/result.out"

if (
    cd "${npz_case_dir}"
    bash ../../integrate/tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a missing NPZ file while generating result.ref" >&2
    exit 1
fi

if [ -e "${npz_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial NPZ result.ref after generation failure" >&2
    exit 1
fi

hsk_case_dir="${test_root}/tests/01_PW/hsk_case"
mkdir -p "${hsk_case_dir}/OUT.autotest"

cat > "${hsk_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
out_hsk 2
gamma_only 1
EOF

cat > "${hsk_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${hsk_case_dir}/sr_nao.csr.ref"
touch "${hsk_case_dir}/OUT.autotest/sr_nao.csr"

(
    cd "${hsk_case_dir}"
    bash -e ../../integrate/tools/catch_properties.sh result.out
)

grep -qx 'CompareH_pass 1' "${hsk_case_dir}/result.out"
grep -qx 'CompareS_pass 1' "${hsk_case_dir}/result.out"

if (
    cd "${hsk_case_dir}"
    bash ../../integrate/tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a failed HSK comparison while generating result.ref" >&2
    exit 1
fi

if [ -e "${hsk_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial HSK result.ref after generation failure" >&2
    exit 1
fi

cat > "${test_root}/fake_abacus" <<'EOF'
#!/usr/bin/env bash
mkdir -p OUT.autotest
cat > OUT.autotest/running_get_s.log <<'INNER_EOF'
Driver atomic_world 0.01
INNER_EOF
touch OUT.autotest/sr_nao.csr
EOF
chmod +x "${test_root}/fake_abacus"

echo test_case > "${test_root}/tests/01_PW/CASES.txt"
autotest_log="${test_root}/autotest.log"
if (
    cd "${test_root}/tests/01_PW"
    OMP_NUM_THREADS=1 bash ../integrate/Autotest.sh \
        -a "${test_root}/fake_abacus" -n 0 -o 1 -j 1 -f CASES.txt -g
) >"${autotest_log}" 2>&1; then
    cat "${autotest_log}" >&2
    echo "Autotest.sh reported success after reference generation failed" >&2
    exit 1
fi

if [ -e "${case_dir}/result.ref" ]; then
    echo "Autotest.sh left a partial result.ref after generation failure" >&2
    exit 1
fi
