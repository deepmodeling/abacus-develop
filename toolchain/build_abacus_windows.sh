#!/bin/bash -e
# Build ABACUS natively on Windows (MSYS2 / MinGW-w64), serial plane-wave only.
#
# Windows counterpart of build_abacus_gnu.sh. Run it from the "MSYS2 MinGW
# 64-bit" shell after ./toolchain_windows.sh has installed the prerequisites.

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
[ -f "$INSTALL_DIR/setup" ] && source "$INSTALL_DIR/setup"
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)

BUILD_DIR=build_abacus_windows
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
LAPACK=${OPENBLAS_ROOT}/lib   # OpenBLAS supplies both BLAS and LAPACK
FFTW3=${FFTW_ROOT}

NUM_JOBS="$(nproc)"
while [[ $# -gt 0 ]]; do
  case $1 in
    -j)
      if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then NUM_JOBS="${2}"; shift 2
      else echo "ERROR: -j requires a number argument"; exit 1; fi ;;
    -j[0-9]*) NUM_JOBS="${1#-j}"; shift ;;
    *) echo "ERROR: Unsupported argument: $1" >&2; echo "Usage: $0 [-j N|-jN]" >&2; exit 1 ;;
  esac
done

# Notes on the non-default options:
#  * ENABLE_MPI/LCAO/ELPA/PEXSI/LIBRI/MLALGO/CUDA = OFF -> Phase 1 serial PW.
#  * BLA_VENDOR=OpenBLAS          -> let CMake's FindBLAS/FindLAPACK pick OpenBLAS.
#  * ENABLE_FLOAT_FFTW=ON         -> make FFT_CPU<float> concrete (vtable) on PE.
#  * COMMIT_INFO=OFF              -> skip the git/sh build-stamp step.
#  * CMAKE_CXX_FLAGS "-include .." -> MSYS2 ships a very new GCC whose libstdc++
#      dropped transitive standard headers; force-include the common ones so the
#      existing sources build unchanged. (Not Windows-specific; tied to GCC>=15.)
cmake -B $BUILD_DIR -G Ninja -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=gcc \
        -DCMAKE_CXX_COMPILER=g++ \
        -DENABLE_MPI=OFF \
        -DENABLE_LCAO=OFF \
        -DUSE_OPENMP=OFF \
        -DUSE_ELPA=OFF \
        -DENABLE_PEXSI=OFF \
        -DENABLE_LIBRI=OFF \
        -DENABLE_MLALGO=OFF \
        -DUSE_CUDA=OFF \
        -DBUILD_TESTING=OFF \
        -DCOMMIT_INFO=OFF \
        -DBLA_VENDOR=OpenBLAS \
        -DENABLE_FLOAT_FFTW=ON \
        -DLAPACK_DIR=$LAPACK \
        -DFFTW3_DIR=$FFTW3 \
        -DCMAKE_PREFIX_PATH=${MINGW_PREFIX:-/mingw64} \
        -DCMAKE_CXX_FLAGS="-include cstdint -include cstring -include algorithm"

cmake --build $BUILD_DIR -j "${NUM_JOBS}"

# Provide a generic `abacus` command, matching the Linux toolchain (which
# symlinks `abacus` -> abacus_<config>). Native Windows symlinks need elevated
# privileges, so instead copy the built binary to abacus.exe; a bare `abacus`
# then resolves to it in the MSYS2 shell (and in cmd/PowerShell). The glob
# matches the configured target (abacus_pw_ser.exe now, abacus_basic_ser.exe in
# later phases) but not the abacus.exe copy itself (no underscore).
built_exe=$(ls "${ABACUS_DIR}/${BUILD_DIR}"/abacus_*.exe 2>/dev/null | head -n 1)
if [ -n "$built_exe" ]; then
    cp -f "$built_exe" "${ABACUS_DIR}/${BUILD_DIR}/abacus.exe"
    echo "Created generic launcher: ${ABACUS_DIR}/${BUILD_DIR}/abacus.exe -> $(basename "$built_exe")"
else
    echo "WARNING: no abacus_*.exe found in ${BUILD_DIR}; 'abacus' command not created."
fi

# generate abacus_env.sh: sourcing it puts the MinGW runtime DLLs (via the
# toolchain setup) and the binary directory on PATH, so `abacus` runs directly.
cat << EOF > "${TOOL}/abacus_env.sh"
#!/bin/bash
[ -f "${INSTALL_DIR}/setup" ] && source "${INSTALL_DIR}/setup"
export PATH="${ABACUS_DIR}/${BUILD_DIR}":\${PATH}
EOF

cat << EOF
========================== usage =========================
Done! Binary: ${ABACUS_DIR}/${BUILD_DIR}/abacus_pw_ser.exe
To run it, exactly like the Linux toolchain in a mingw bash:
    bash
    source ${TOOL}/abacus_env.sh
    abacus                 # -> abacus.exe (a copy of abacus_pw_ser.exe)
==========================================================
EOF
