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
    bash -e ../../integrate/validation_tools/catch_properties.sh result.out
)

grep -qx 'CompareBand_pass 1' "${case_dir}/result.out"

if (
    cd "${case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
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
    bash -e ../../integrate/validation_tools/catch_properties.sh result.out
)

grep -qx 'OutputSRNPZ_pass 1' "${npz_case_dir}/result.out"
grep -qx 'OutputHRNPZ_pass 1' "${npz_case_dir}/result.out"

if (
    cd "${npz_case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
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
    bash -e ../../integrate/validation_tools/catch_properties.sh result.out
)

grep -qx 'CompareH_pass 1' "${hsk_case_dir}/result.out"
grep -qx 'CompareS_pass 1' "${hsk_case_dir}/result.out"

if (
    cd "${hsk_case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a failed HSK comparison while generating result.ref" >&2
    exit 1
fi

if [ -e "${hsk_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial HSK result.ref after generation failure" >&2
    exit 1
fi

# Keep upstream's renamed VXC files while retaining failure propagation in
# both Gamma-only and multi-k collectors. No legacy vxc_out.* files exist.
for gamma in 0 1; do
    vxc_case_dir="${test_root}/tests/01_PW/vxc_case_${gamma}"
    mkdir -p "${vxc_case_dir}/OUT.autotest"
    printf 'INPUT_PARAMETERS\ncalculation get_s\nout_mat_xc 1\ngamma_only %s\n' \
        "${gamma}" > "${vxc_case_dir}/INPUT"
    printf 'Driver atomic_world 0.01\n' > "${vxc_case_dir}/OUT.autotest/running_get_s.log"
    touch "${vxc_case_dir}/sr_nao.csr.ref" "${vxc_case_dir}/OUT.autotest/sr_nao.csr"
    if [ "${gamma}" -eq 1 ]; then
        vxc_name=vxc_nao.txt
    else
        vxc_name=vxck2_nao.txt
    fi
    printf 'value 1\n' > "${vxc_case_dir}/${vxc_name}.ref"
    printf 'value 1\n' > "${vxc_case_dir}/OUT.autotest/${vxc_name}"
    printf 'value 1\n' > "${vxc_case_dir}/vxc.ref"
    printf 'value 1\n' > "${vxc_case_dir}/OUT.autotest/vxc.txt"
    (
        cd "${vxc_case_dir}"
        bash -e ../../integrate/validation_tools/catch_properties.sh result.ref
    )
    grep -qx 'CompareVXC_pass 0' "${vxc_case_dir}/result.ref"
    grep -qx 'CompareOrbXC_pass 0' "${vxc_case_dir}/result.ref"

    for failure in different missing; do
        if [ "${failure}" = different ]; then
            printf 'value 2\n' > "${vxc_case_dir}/OUT.autotest/vxc.txt"
        else
            rm "${vxc_case_dir}/OUT.autotest/vxc.txt"
        fi
        (
            cd "${vxc_case_dir}"
            bash -e ../../integrate/validation_tools/catch_properties.sh result.out
        )
        grep -qx 'CompareVXC_pass 0' "${vxc_case_dir}/result.out"
        grep -qx 'CompareOrbXC_pass 1' "${vxc_case_dir}/result.out"
        if (
            cd "${vxc_case_dir}"
            bash -e ../../integrate/validation_tools/catch_properties.sh result.ref
        ); then
            echo "catch_properties.sh accepted ${failure} VXC output in result.ref" >&2
            exit 1
        fi
        if [ -e "${vxc_case_dir}/result.ref" ]; then
            echo "catch_properties.sh left a partial VXC result.ref" >&2
            exit 1
        fi
    done
done

# Upstream's Wannier validators accept the band gauge without rewriting the
# output files; all three statuses must still reject invalid references.
wannier_case_dir="${test_root}/tests/01_PW/wannier_case"
mkdir -p "${wannier_case_dir}/OUT.autotest"
printf 'INPUT_PARAMETERS\ncalculation get_s\ntowannier90 1\n' > "${wannier_case_dir}/INPUT"
printf 'Driver atomic_world 0.01\n' > "${wannier_case_dir}/OUT.autotest/running_get_s.log"
touch "${wannier_case_dir}/sr_nao.csr.ref" "${wannier_case_dir}/OUT.autotest/sr_nao.csr"
printf '2 1 1\n1 1 1 1 0\n2 1 1 0 0\n' > "${wannier_case_dir}/diamond.amn"
printf 'Created on test date\n2 1 1\n1 1 1 0 0\n2 1 1 0 1\n' > "${wannier_case_dir}/OUT.autotest/diamond.amn"
printf '1 1 1\n1 1 0 0 0\n1 0\n' > "${wannier_case_dir}/diamond.mmn"
printf 'Created on test date\n1 1 1\n1 1 0 0 0\n0 1\n' > "${wannier_case_dir}/OUT.autotest/diamond.mmn"
printf '1 1 0.5\n' > "${wannier_case_dir}/diamond.eig"
cp "${wannier_case_dir}/diamond.eig" "${wannier_case_dir}/OUT.autotest/diamond.eig"
(
    cd "${wannier_case_dir}"
    bash -e ../../integrate/validation_tools/catch_properties.sh result.ref
)
for property in AMN MMN EIG; do
    grep -qx "Compare${property}_pass 0" "${wannier_case_dir}/result.ref"
done
grep -qx 'Created on test date' "${wannier_case_dir}/OUT.autotest/diamond.amn"
grep -qx 'Created on test date' "${wannier_case_dir}/OUT.autotest/diamond.mmn"
for extension in amn mmn eig; do
    property=${extension^^}
    output_file="${wannier_case_dir}/OUT.autotest/diamond.${extension}"
    cp "${output_file}" "${test_root}/valid.${extension}"
    for failure in malformed missing; do
        if [ "${failure}" = malformed ]; then
            printf 'invalid Wannier output\n' > "${output_file}"
        else
            rm "${output_file}"
        fi
        (
            cd "${wannier_case_dir}"
            bash -e ../../integrate/validation_tools/catch_properties.sh result.out
        )
        grep -qx "Compare${property}_pass 1" "${wannier_case_dir}/result.out"
        if (
            cd "${wannier_case_dir}"
            bash ../../integrate/validation_tools/catch_properties.sh result.ref
        ); then
            echo "catch_properties.sh accepted ${failure} ${property} output in result.ref" >&2
            exit 1
        fi
        test ! -e "${wannier_case_dir}/result.ref"
    done
    cp "${test_root}/valid.${extension}" "${output_file}"
done

# Both upstream LR routes keep failed checks in result.out, but reference
# generation must fail and remove even the early-return tracking result.
for lr_mode in force tracking; do
    lr_case_dir="${test_root}/tests/01_PW/lr_${lr_mode}"
    mkdir -p "${lr_case_dir}/OUT.autotest"
    printf 'INPUT_PARAMETERS\ncalculation get_s\nesolver_type lr\ncal_force 1\n' > "${lr_case_dir}/INPUT"
    printf 'Driver atomic_world 0.01\n' > "${lr_case_dir}/OUT.autotest/running_get_s.log"
    touch "${lr_case_dir}/sr_nao.csr.ref" "${lr_case_dir}/OUT.autotest/sr_nao.csr"
    if [ "${lr_mode}" = force ]; then
        lr_property=CompareLRForce_pass
        printf 'singlet 0 0 1 -2 3\n' > "${lr_case_dir}/lr_force.ref"
        cat >> "${lr_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Start to calculate excited-state force of singlet
Forces (-gradients) of each excited state: (eV/Angstrom)
 state atom x y z
 0 0 force 1 -2 3
DONE
EOF
    else
        lr_property=CompareStateTracking_pass
        printf 'out_band 1\n' >> "${lr_case_dir}/INPUT"
        printf '{"steps": 1}\n' > "${lr_case_dir}/state_tracking.ref"
        cat >> "${lr_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
EXCITED-STATE RELAX step 0: E_gs = -10 eV, Omega = 1 eV, E_exc = -9 eV | |F_gs|max = 0.1, |F_Omega|max = 0.2 eV/Angstrom
EOF
    fi
    (
        cd "${lr_case_dir}"
        bash -e ../../integrate/validation_tools/catch_properties.sh result.ref
    )
    grep -qx "${lr_property} 0" "${lr_case_dir}/result.ref"
    if [ "${lr_mode}" = tracking ]; then
        test "$(wc -l < "${lr_case_dir}/result.ref")" -eq 1
    fi
    printf 'Driver atomic_world 0.01\n' > "${lr_case_dir}/OUT.autotest/running_get_s.log"
    (
        cd "${lr_case_dir}"
        bash -e ../../integrate/validation_tools/catch_properties.sh result.out
    )
    grep -qx "${lr_property} 1" "${lr_case_dir}/result.out"
    if (
        cd "${lr_case_dir}"
        bash ../../integrate/validation_tools/catch_properties.sh result.ref
    ); then
        echo "catch_properties.sh accepted a failed LR ${lr_mode} check in result.ref" >&2
        exit 1
    fi
    test ! -e "${lr_case_dir}/result.ref"
done

deepks_case_dir="${test_root}/tests/01_PW/deepks_case"
mkdir -p "${deepks_case_dir}/OUT.autotest"

cat > "${deepks_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
deepks_out_labels 1
deepks_v_delta -1
EOF

cat > "${deepks_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${deepks_case_dir}/sr_nao.csr.ref"
touch "${deepks_case_dir}/OUT.autotest/sr_nao.csr"
printf 'value 1\n' > "${deepks_case_dir}/deepks_hrtot.csr.ref"
printf 'value 2\n' > "${deepks_case_dir}/OUT.autotest/deepks_hrtot.csr"
printf 'value 1\n' > "${deepks_case_dir}/deepks_hrdelta.csr.ref"
printf 'value 1\n' > "${deepks_case_dir}/OUT.autotest/deepks_hrdelta.csr"

(
    cd "${deepks_case_dir}"
    bash -e ../../integrate/validation_tools/catch_properties.sh result.out
)

grep -qx 'deepks_hr_label_pass 1' "${deepks_case_dir}/result.out"
grep -qx 'deepks_vdelta_r_pass 0' "${deepks_case_dir}/result.out"

if (
    cd "${deepks_case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a failed DeePKS comparison while generating result.ref" >&2
    exit 1
fi

if [ -e "${deepks_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial DeePKS result.ref after generation failure" >&2
    exit 1
fi

deepks_missing_case_dir="${test_root}/tests/01_PW/deepks_missing_case"
mkdir -p "${deepks_missing_case_dir}/OUT.autotest"

cat > "${deepks_missing_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
deepks_scf 1
EOF

cat > "${deepks_missing_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${deepks_missing_case_dir}/sr_nao.csr.ref"
touch "${deepks_missing_case_dir}/OUT.autotest/sr_nao.csr"

(
    cd "${deepks_missing_case_dir}"
    bash -e ../../integrate/validation_tools/catch_properties.sh result.out
)

grep -qx 'deepks_desc 0.0' "${deepks_missing_case_dir}/result.out"

if (
    cd "${deepks_missing_case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a missing DeePKS descriptor while generating result.ref" >&2
    exit 1
fi

if [ -e "${deepks_missing_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial result.ref after a missing DeePKS descriptor" >&2
    exit 1
fi

ml_case_dir="${test_root}/tests/01_PW/ml_case"
mkdir -p "${ml_case_dir}/OUT.autotest/MLKEDF_Descriptors"

cat > "${ml_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
EOF

cat > "${ml_case_dir}/OUT.autotest/running_get_s.log" <<'EOF'
Driver atomic_world 0.01
EOF

touch "${ml_case_dir}/sr_nao.csr.ref"
touch "${ml_case_dir}/OUT.autotest/sr_nao.csr"
printf 'not a NumPy file\n' > "${ml_case_dir}/OUT.autotest/MLKEDF_Descriptors/broken.npy"

if (
    cd "${ml_case_dir}"
    bash ../../integrate/validation_tools/catch_properties.sh result.ref
); then
    echo "catch_properties.sh accepted a failed ML descriptor collector" >&2
    exit 1
fi

if [ -e "${ml_case_dir}/result.ref" ]; then
    echo "catch_properties.sh left a partial result.ref after an ML collector failure" >&2
    exit 1
fi

cat > "${test_root}/fake_abacus" <<'EOF'
#!/usr/bin/env bash
mkdir -p OUT.autotest
cat > OUT.autotest/running_get_s.log <<'INNER_EOF'
Driver atomic_world 0.01
INNER_EOF
cat > OUT.autotest/running_scf.log <<'INNER_EOF'
Driver atomic_world 0.01
INNER_EOF
cat > OUT.autotest/running_get_wf.log <<'INNER_EOF'
Driver atomic_world 0.01
INNER_EOF
touch OUT.autotest/sr_nao.csr
EOF
chmod +x "${test_root}/fake_abacus"

mkdir -p "${test_root}/fake_bin"
cat > "${test_root}/fake_bin/mpirun" <<'EOF'
#!/usr/bin/env bash
if [ "$1" = "-np" ]; then
    shift 2
fi
exec "$@"
EOF
chmod +x "${test_root}/fake_bin/mpirun"

# A collector can fail after writing a valid prefix of result.out.  run_check.sh
# must propagate that failure before check_out can accept the prefix. Exercise
# both the ordinary path and Single_job.sh's debug path, which uses the copied
# validation_tools directory.
for run_check_mode in normal debug; do
    run_check_case_dir="${test_root}/tests/01_PW/run_check_${run_check_mode}"
    mkdir -p "${run_check_case_dir}/OUT.autotest/MLKEDF_Descriptors"
    printf 'not a NumPy file\n' > "${run_check_case_dir}/OUT.autotest/MLKEDF_Descriptors/broken.npy"
    cat > "${run_check_case_dir}/INPUT" <<'EOF'
INPUT_PARAMETERS
calculation get_s
EOF
    printf 'Driver atomic_world 0.01\n' > "${run_check_case_dir}/OUT.autotest/running_get_s.log"
    touch "${run_check_case_dir}/sr_nao.csr.ref" "${run_check_case_dir}/OUT.autotest/sr_nao.csr"
    printf 'CompareS_pass 0\n' > "${run_check_case_dir}/result.ref"

    if [ "${run_check_mode}" = debug ]; then
        cp "${integrate_dir}/validation_tools/run_check.sh" "${run_check_case_dir}/run_check.sh"
        cp -r "${integrate_dir}/validation_tools" "${run_check_case_dir}/validation_tools"
        run_check_command=(bash ./run_check.sh debug)
    else
        run_check_command=(bash ../../integrate/validation_tools/run_check.sh)
    fi

    run_check_log="${test_root}/run-check-${run_check_mode}.log"
    if (
        cd "${run_check_case_dir}"
        PATH="${test_root}/fake_bin:${PATH}" \
            ABACUS_EXE="${test_root}/fake_abacus" \
            "${run_check_command[@]}"
    ) >"${run_check_log}" 2>&1; then
        cat "${run_check_log}" >&2
        echo "run_check.sh ${run_check_mode} path reported success after collector failure" >&2
        exit 1
    fi
    grep -q 'Failed to collect validation results.' "${run_check_log}"
done

run_check_log="${test_root}/run-check.log"
if (
    cd "${deepks_missing_case_dir}"
    PATH="${test_root}/fake_bin:${PATH}" \
        ABACUS_EXE="${test_root}/fake_abacus" \
        bash ../../integrate/validation_tools/run_check.sh generate
) >"${run_check_log}" 2>&1; then
    cat "${run_check_log}" >&2
    echo "run_check.sh reported success after reference generation failed" >&2
    exit 1
fi

if [ -e "${deepks_missing_case_dir}/result.ref" ]; then
    echo "run_check.sh left a partial result.ref after generation failure" >&2
    exit 1
fi

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

# Upstream's dedicated-validator path must still honor a custom MPI launcher,
# its flags, the OFML single-rank exception, and commented case-list entries.
launcher_case_dir="${test_root}/tests/01_PW/launcher_case"
mkdir -p "${launcher_case_dir}"
printf 'INPUT_PARAMETERS\ncalculation get_s\nout_band 1\n' > "${launcher_case_dir}/INPUT"
printf '\n  # absent_case\nlauncher_case\n' > "${test_root}/tests/01_PW/LAUNCHER_CASES.txt"
cat > "${test_root}/fake_bin/custom_mpi" <<'EOF'
#!/usr/bin/env bash
set -eu
[ "$#" -eq 6 ]
[ "$1" = "--ranks" ]
[ "$2" = "$EXPECTED_RANKS" ]
[ "$3" = "--test-pre-one" ]
[ "$4" = "--test-pre-two" ]
[ "$6" = "--test-post" ]
touch "$LAUNCHER_MARKER"
exec "$5"
EOF
chmod +x "${test_root}/fake_bin/custom_mpi"
for expected_ranks in 4 1; do
    if [ "${expected_ranks}" -eq 1 ]; then
        printf 'of_ml_gene_data 1\n' >> "${launcher_case_dir}/INPUT"
    fi
    launcher_log="${test_root}/launcher-${expected_ranks}.log"
    (
        cd "${test_root}/tests/01_PW"
        EXPECTED_RANKS="${expected_ranks}" \
        LAUNCHER_MARKER="${test_root}/launcher-${expected_ranks}.used" \
        ABACUS_MPIEXEC="${test_root}/fake_bin/custom_mpi" \
        ABACUS_MPIEXEC_NUMPROC_FLAG=--ranks \
        ABACUS_MPIEXEC_PREFLAGS='--test-pre-one --test-pre-two' \
        ABACUS_MPIEXEC_POSTFLAGS=--test-post \
            OMP_NUM_THREADS=1 bash ../integrate/Autotest.sh \
            -a "${test_root}/fake_abacus" -n 4 -o 1 -j 1 -f LAUNCHER_CASES.txt -p
    ) > "${launcher_log}" 2>&1
    test -f "${test_root}/launcher-${expected_ranks}.used"
    grep -q 'Property collectors skipped' "${launcher_log}"
    if [ -e "${launcher_case_dir}/result.out" ] || [ -e "${launcher_case_dir}/result.ref" ]; then
        echo "Autotest.sh ran property collectors despite -p" >&2
        exit 1
    fi
done
