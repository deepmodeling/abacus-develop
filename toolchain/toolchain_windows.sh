#!/bin/bash -e
# Toolchain setup for a NATIVE Windows build of ABACUS via MSYS2 / MinGW-w64.
#
# This is the Windows counterpart of toolchain_gnu.sh / toolchain_intel.sh.
# On Linux those scripts build the dependency stack from source; on Windows the
# MinGW-w64 dependencies are provided by the MSYS2 distribution, so here we just
# install them with pacman and record their location for build_abacus_windows.sh.
#
# Scope: PW + LCAO, serial and MPI (MS-MPI + ScaLAPACK). ELPA / PEXSI / hybrid
# functionals (LibRI) / DeePKS / LibXC / GPU are intentionally omitted because
# they have no reliable native-Windows build yet; they remain ordinary ABACUS
# feature switches for the future.
#
# Usage: open the "MSYS2 MinGW 64-bit" shell and run:
#     ./toolchain_windows.sh
# then:
#     ./build_abacus_windows.sh

if ! command -v pacman >/dev/null 2>&1; then
    echo "ERROR: pacman not found. Run this inside an MSYS2 shell (https://www.msys2.org)."
    exit 1
fi

echo "[*] Installing MinGW-w64 build prerequisites via pacman ..."
pacman -S --needed --noconfirm \
    mingw-w64-x86_64-gcc \
    mingw-w64-x86_64-gcc-fortran \
    mingw-w64-x86_64-cmake \
    mingw-w64-x86_64-ninja \
    mingw-w64-x86_64-openblas \
    mingw-w64-x86_64-fftw \
    mingw-w64-x86_64-cereal \
    mingw-w64-x86_64-msmpi \
    mingw-w64-x86_64-scalapack

# Notes:
#  * cereal    : header-only serialization, required by the LCAO build.
#  * msmpi     : MS-MPI headers + import lib for the MPI build. The MS-MPI
#                *runtime* (msmpi.dll, mpiexec) is a separate Microsoft
#                redistributable that must be installed system-wide to run
#                parallel jobs: https://www.microsoft.com/download/details.aspx?id=105289
#  * scalapack : distributed eigensolver used by the LCAO MPI build (no ELPA).

# 'bc' (a base MSYS tool, not a MinGW package) is used by the integration-test
# harness tests/integrate/tools/catch_properties.sh; install it so the existing
# serial test flow (Autotest.sh -n 0) works out of the box.
pacman -S --needed --noconfirm bc

# MinGW-w64 installs everything under the /mingw64 prefix. Record it in a setup
# file with the same variable names the build_abacus_*.sh scripts expect, so the
# build step is uniform with the Linux toolchain.
TOOL=$(cd "$(dirname "$0")" && pwd)
INSTALL_DIR="$TOOL/install"
mkdir -p "$INSTALL_DIR"
cat > "$INSTALL_DIR/setup" <<'EOF'
# Native Windows (MSYS2/MinGW-w64) prerequisites live under /mingw64.
export MINGW_PREFIX="${MINGW_PREFIX:-/mingw64}"
export OPENBLAS_ROOT="$MINGW_PREFIX"   # OpenBLAS provides BLAS *and* LAPACK
export FFTW_ROOT="$MINGW_PREFIX"
export PATH="$MINGW_PREFIX/bin:$PATH"
EOF

cat <<EOF
========================== done ==========================
MinGW-w64 prerequisites installed; environment recorded in:
    $INSTALL_DIR/setup
Next:
    ./build_abacus_windows.sh
==========================================================
EOF
