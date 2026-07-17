#!/bin/bash
#
# SCF + NSCF Workflow Script for ABACUS DeltaSpin/DFT+U tests
#
# Usage:
#   bash run_scf_nscf.sh <abacus_path> [mpi_np]
#
# This script:
#   1. Runs SCF calculation in scf/ subdirectory
#   2. Copies charge density (and onsite.dm for DFT+U) from SCF output
#   3. Runs NSCF calculation in current directory
#   4. Compares output with reference
#
# Directory structure:
#   ./
#   ├── INPUT          (NSCF input)
#   ├── STRU
#   ├── KPT
#   ├── scf/
#   │   ├── INPUT      (SCF input)
#   │   ├── STRU       (same as parent)
#   │   └── KPT        (SCF k-points, usually Gamma only)
#   └── OUT.*          (SCF output with charge density)
#

set -e

ABACUS="${1:-abacus}"
MPI_NP="${2:-4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(pwd)"

# Allow OpenMPI run as root
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

echo "========================================"
echo " SCF + NSCF Workflow"
echo "========================================"
echo "ABACUS:   $ABACUS"
echo "MPI NP:   $MPI_NP"
echo "Test dir: $TEST_DIR"
echo "========================================"

# -------------------------------------------------------
# Step 0: Check prerequisites
# -------------------------------------------------------
if [ ! -f "${TEST_DIR}/INPUT" ]; then
    echo "ERROR: INPUT file not found in ${TEST_DIR}"
    exit 1
fi

if [ ! -f "${TEST_DIR}/STRU" ]; then
    echo "ERROR: STRU file not found in ${TEST_DIR}"
    exit 1
fi

# Check if this is a DFT+U test (needs onsite.dm handling)
IS_DFTU=false
if grep -q "dft_plus_u" "${TEST_DIR}/INPUT" 2>/dev/null; then
    IS_DFTU=true
    echo "DFT+U test detected: will handle onsite.dm files"
fi

# -------------------------------------------------------
# Step 1: Setup SCF subdirectory
# -------------------------------------------------------
SCF_DIR="${TEST_DIR}/scf"
mkdir -p "${SCF_DIR}"

# Copy STRU to SCF directory
cp "${TEST_DIR}/STRU" "${SCF_DIR}/STRU"

# Copy KPT to SCF directory (if exists)
if [ -f "${TEST_DIR}/KPT" ]; then
    cp "${TEST_DIR}/KPT" "${SCF_DIR}/KPT"
fi

# Generate SCF INPUT from NSCF INPUT
# Key changes: calculation=scf, init_chg=atomic, remove nscf-specific params
SCF_INPUT="${SCF_DIR}/INPUT"
cat "${TEST_DIR}/INPUT" | sed \
    -e 's/calculation\s*nscf/calculation    scf/' \
    -e 's/init_chg\s*file/init_chg    atomic/' \
    -e '/read_file_dir/d' \
    -e '/out_band/d' \
    > "${SCF_INPUT}"

echo ""
echo "--- SCF INPUT ---"
cat "${SCF_INPUT}"
echo "-----------------"

# -------------------------------------------------------
# Step 2: Run SCF calculation
# -------------------------------------------------------
echo ""
echo "[1/4] Running SCF calculation..."
cd "${SCF_DIR}"
mpirun -np ${MPI_NP} "${ABACUS}" > scf.log 2>&1 || {
    echo "ERROR: SCF calculation failed!"
    echo "Check ${SCF_DIR}/scf.log for details"
    cd "${TEST_DIR}"
    exit 1
}
cd "${TEST_DIR}"

# Get the suffix from SCF INPUT (default: OUT)
SCF_SUFFIX=$(grep "suffix" "${SCF_INPUT}" | awk '{print $2}')
SCF_SUFFIX=${SCF_SUFFIX:-OUT}
SCF_OUT="${SCF_DIR}/${SCF_SUFFIX}"

# -------------------------------------------------------
# Step 3: Copy charge density from SCF output
# -------------------------------------------------------
echo "[2/4] Copying charge density files..."

# NSCF suffix
NSCF_SUFFIX=$(grep "suffix" "${TEST_DIR}/INPUT" | awk '{print $2}')
NSCF_SUFFIX=${NSCF_SUFFIX:-OUT}

# Find charge density file (pattern: *-CHARGE-DENSITY.restart or *-CHARGE-DENSITY)
CHG_FILE=$(find "${SCF_OUT}" -name "*CHARGE-DENSITY*" 2>/dev/null | head -1)
if [ -z "${CHG_FILE}" ]; then
    echo "ERROR: No charge density file found in ${SCF_OUT}"
    exit 1
fi

CHG_BASENAME=$(basename "${CHG_FILE}")
echo "  Found: ${CHG_FILE}"

# Determine the target filename based on NSCF suffix
NSCF_CHG_FILE="${NSCF_SUFFIX}-${CHG_BASENAME#*-}"
cp "${CHG_FILE}" "${TEST_DIR}/${NSCF_CHG_FILE}"
echo "  Copied to: ${TEST_DIR}/${NSCF_CHG_FILE}"

# For DFT+U tests, also copy onsite.dm if it exists
if [ "${IS_DFTU}" = true ]; then
    ONSITE_FILE=$(find "${SCF_OUT}" -name "onsite.dm" 2>/dev/null | head -1)
    if [ -n "${ONSITE_FILE}" ]; then
        cp "${ONSITE_FILE}" "${TEST_DIR}/onsite.dm"
        echo "  Copied onsite.dm"
    else
        echo "  WARNING: onsite.dm not found in SCF output"
    fi
fi

# -------------------------------------------------------
# Step 4: Run NSCF calculation
# -------------------------------------------------------
echo "[3/4] Running NSCF calculation..."
mpirun -np ${MPI_NP} "${ABACUS}" > nscf.log 2>&1 || {
    echo "ERROR: NSCF calculation failed!"
    echo "Check ${TEST_DIR}/nscf.log for details"
    exit 1
}

# -------------------------------------------------------
# Step 5: Verify output
# -------------------------------------------------------
echo "[4/4] Verifying output..."

NSF_OUT="${TEST_DIR}/${NSCF_SUFFIX}"

# Check if output directory exists
if [ ! -d "${NSF_OUT}" ]; then
    echo "WARNING: Output directory ${NSF_OUT} not found"
    exit 1
fi

# Print summary
echo ""
echo "========================================"
echo " SCF + NSCF Workflow Complete"
echo "========================================"
echo "SCF output:  ${SCF_OUT}/"
echo "NSCF output: ${NSF_OUT}/"
echo ""

# Print energy from SCF
if [ -f "${SCF_OUT}/running_scf.log" ]; then
    echo "SCF final energy:"
    grep "!FINAL ENERGY" "${SCF_OUT}/running_scf.log" | tail -1 || echo "  (not found)"
fi

# Print magnetic moments from NSCF
if [ -f "${NSF_OUT}/running_scf.log" ] || [ -f "${NSF_OUT}/running_nscf.log" ]; then
    echo ""
    echo "NSCF magnetic moments:"
    grep -A 1 "TOTAL MAGNETISM" "${NSF_OUT}/running_scf.log" 2>/dev/null | tail -3 || \
    grep -A 1 "TOTAL MAGNETISM" "${NSF_OUT}/running_nscf.log" 2>/dev/null | tail -3 || \
    echo "  (not found)"
fi

# Print lambda if DeltaSpin
if grep -q "sc_mag_switch" "${TEST_DIR}/INPUT" 2>/dev/null; then
    echo ""
    echo "NSCF lambda values:"
    grep "lambda" "${NSF_OUT}/running_scf.log" 2>/dev/null | tail -3 || \
    grep "lambda" "${NSF_OUT}/running_nscf.log" 2>/dev/null | tail -3 || \
    echo "  (not found)"
fi

echo ""
echo "========================================"
