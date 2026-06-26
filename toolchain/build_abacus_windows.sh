#!/bin/bash -e
# Build ABACUS natively on Windows (MSYS2 / MinGW-w64).
#
# Windows counterpart of build_abacus_gnu.sh. Run it from the "MSYS2 MinGW
# 64-bit" shell after ./toolchain_windows.sh has installed the prerequisites.
#
# By default it builds the most capable supported configuration: MPI + LCAO
# (plane-wave and numerical-atomic-orbital bases) with OpenBLAS + FFTW +
# ScaLAPACK. ELPA / PEXSI / hybrid functionals (LibRI) / DeePKS / GPU are not
# available on Windows yet and stay OFF.
#
# Override the configuration from the environment, e.g.:
#   ENABLE_MPI=OFF ./build_abacus_windows.sh     # serial
#   ENABLE_LCAO=OFF ./build_abacus_windows.sh    # plane-wave only
#   ENABLE_MPI=OFF ENABLE_LCAO=OFF ./build_abacus_windows.sh  # serial PW (Phase 1)
ENABLE_MPI=${ENABLE_MPI:-ON}
ENABLE_LCAO=${ENABLE_LCAO:-ON}

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
[ -f "$INSTALL_DIR/setup" ] && source "$INSTALL_DIR/setup"
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)
MINGW_PREFIX=${MINGW_PREFIX:-/mingw64}

BUILD_DIR=build_abacus_windows
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
LAPACK=${OPENBLAS_ROOT:-$MINGW_PREFIX}/lib   # OpenBLAS supplies both BLAS and LAPACK
FFTW3=${FFTW_ROOT:-$MINGW_PREFIX}

NUM_JOBS="$(nproc)"
# Cap the *default* parallelism by available RAM. Several heavy -O3 template
# TUs (e.g. source_cell/module_symmetry/symmetry.cpp, read_pp_upf201.cpp) can
# each peak around 3 GB in cc1plus, and ninja tends to schedule them together;
# on a many-core box -j nproc then exhausts memory and the build dies with
# "cc1plus.exe: out of memory" (seen even on a 31 GB / 20-core machine at
# -j 20). Budget ~3 GB per job. An explicit -j below always overrides this.
if [ -r /proc/meminfo ]; then
  mem_gb=$(awk '/^MemTotal:/ {printf "%d", $2/1024/1024}' /proc/meminfo)
  if [ -n "$mem_gb" ] && [ "$mem_gb" -ge 1 ]; then
    mem_jobs=$(( mem_gb / 3 )); [ "$mem_jobs" -lt 1 ] && mem_jobs=1
    [ "$mem_jobs" -lt "$NUM_JOBS" ] && NUM_JOBS=$mem_jobs
  fi
fi
while [[ $# -gt 0 ]]; do
  case $1 in
    -j)
      if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then NUM_JOBS="${2}"; shift 2
      else echo "ERROR: -j requires a number argument"; exit 1; fi ;;
    -j[0-9]*) NUM_JOBS="${1#-j}"; shift ;;
    *) echo "ERROR: Unsupported argument: $1" >&2; echo "Usage: $0 [-j N|-jN]" >&2; exit 1 ;;
  esac
done
echo "Building with -j ${NUM_JOBS} (override with -j N; lower it if cc1plus runs out of memory)."

# MPI on Windows is MS-MPI (mingw-w64-x86_64-msmpi). Point FindMPI at it.
MPI_ARGS=()
if [ "$ENABLE_MPI" = "ON" ]; then
  MPI_ARGS=(-DMPI_CXX_INCLUDE_PATH=$MINGW_PREFIX/include
            -DMPI_CXX_LIBRARIES=$MINGW_PREFIX/lib/libmsmpi.dll.a)
fi

# Notes on the non-default options:
#  * USE_ELPA/PEXSI/LIBRI/MLALGO/CUDA = OFF -> not available on Windows yet.
#    When ENABLE_MPI=ON the LCAO solver is ScaLAPACK (found automatically);
#    when serial it is LAPACK (DiagoLapack).
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
        -DENABLE_MPI=$ENABLE_MPI \
        -DENABLE_LCAO=$ENABLE_LCAO \
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
        -DCMAKE_PREFIX_PATH=$MINGW_PREFIX \
        "${MPI_ARGS[@]}" \
        -DCMAKE_CXX_FLAGS="-include cstdint -include cstring -include algorithm"

cmake --build $BUILD_DIR -j "${NUM_JOBS}"

# Provide a generic `abacus` command, matching the Linux toolchain (which
# symlinks `abacus` -> abacus_<config>). Native Windows symlinks need elevated
# privileges, so instead copy the built binary to abacus.exe; a bare `abacus`
# then resolves to it in the MSYS2 shell (and in cmd/PowerShell). The glob
# matches the configured target (abacus_basic_para.exe, abacus_pw_ser.exe, ...)
# but not the abacus.exe copy itself (no underscore).
built_exe=$(ls "${ABACUS_DIR}/${BUILD_DIR}"/abacus_*.exe 2>/dev/null | head -n 1)
if [ -n "$built_exe" ]; then
    cp -f "$built_exe" "${ABACUS_DIR}/${BUILD_DIR}/abacus.exe"
    echo "Created generic launcher: ${ABACUS_DIR}/${BUILD_DIR}/abacus.exe -> $(basename "$built_exe")"
else
    echo "WARNING: no abacus_*.exe found in ${BUILD_DIR}; 'abacus' command not created."
fi

# Bundle the dependent MinGW / OpenBLAS / FFTW / ScaLAPACK runtime DLLs next to
# the binary. Windows searches the *application directory* before PATH, so this
# makes abacus.exe self-contained and, crucially, lets it find its DLLs even
# when launched by a process that does not propagate PATH to its children --
# which is exactly what MS-MPI's mpiexec does when the test harness redirects
# stdout to a file ("error while loading shared libraries"). System DLLs
# (msmpi.dll in System32, kernel32, ...) resolve on their own and are skipped.
if [ -n "$built_exe" ]; then
    echo "Bundling dependent DLLs into ${BUILD_DIR}/ ..."
    ldd "${ABACUS_DIR}/${BUILD_DIR}/abacus.exe" 2>/dev/null \
        | awk -v p="$MINGW_PREFIX" '$3 ~ p {print $3}' | sort -u \
        | while read -r dll; do cp -f "$dll" "${ABACUS_DIR}/${BUILD_DIR}/"; done
fi

# When MPI is on, drop an `mpirun` shim next to the binary so the shared test
# harness (which invokes `mpirun -np N`) drives MS-MPI unchanged. MS-MPI ships
# only `mpiexec`; the shim forwards to it and pins the (OpenMP-threaded) BLAS to
# one thread per rank -- otherwise each rank's multithreaded OpenBLAS
# oversubscribes the cores and its buffer allocator fails under several ranks.
if [ "$ENABLE_MPI" = "ON" ]; then
    cat << 'SHIM' > "${ABACUS_DIR}/${BUILD_DIR}/mpirun"
#!/bin/bash
# mpirun -> mpiexec shim for native Windows (MS-MPI). See build_abacus_windows.sh.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
exec mpiexec "$@"
SHIM
    chmod +x "${ABACUS_DIR}/${BUILD_DIR}/mpirun"
    echo "Created mpirun->mpiexec shim: ${ABACUS_DIR}/${BUILD_DIR}/mpirun"
fi

# generate abacus_env.sh: sourcing it puts the MinGW runtime DLLs (via the
# toolchain setup) and the binary directory on PATH, so `abacus` runs directly.
# MSYS2's OpenBLAS is OpenMP-threaded, so OMP_NUM_THREADS (not the often-cited
# OPENBLAS_NUM_THREADS) is what actually caps its threads; pin it to 1 so that
# `mpiexec -n N abacus` doesn't oversubscribe and trip OpenBLAS's buffer
# allocator. (Both are set; OPENBLAS_NUM_THREADS alone has no effect here.)
cat << EOF > "${TOOL}/abacus_env.sh"
#!/bin/bash
[ -f "${INSTALL_DIR}/setup" ] && source "${INSTALL_DIR}/setup"
export PATH="${ABACUS_DIR}/${BUILD_DIR}":\${PATH}
# MS-MPI's mpiexec lives in its own Bin dir (MSMPI_BIN), which the MinGW PATH
# does not inherit; add it so \`mpiexec\` and the mpirun shim resolve.
[ -n "\$MSMPI_BIN" ] && export PATH="\$(cygpath -u "\$MSMPI_BIN")":\${PATH}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
EOF

cat << EOF
========================== usage =========================
Done! Binary: $(basename "$built_exe") in ${ABACUS_DIR}/${BUILD_DIR}/
Run it from a MinGW bash shell:
    bash
    source ${TOOL}/abacus_env.sh
    abacus                                   # serial run
    mpiexec -n 4 abacus                      # parallel run (MS-MPI)

Run the standard test suite (the mpirun->mpiexec shim makes the existing
harness work unchanged):
    cd ${ABACUS_DIR}/tests/01_PW
    bash ../integrate/Autotest.sh -a abacus          # MPI (default np=4)
    bash ../integrate/Autotest.sh -a abacus -n 0     # serial (no launcher)
==========================================================
EOF
