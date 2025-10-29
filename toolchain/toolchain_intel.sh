#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# Users can easily modify these parameters to customize the build

# Before running this script, ensure you have loaded your system packages
# module load mkl mpi compiler

# Compiler Configuration
TOOLCHAIN_COMPILER="intel"
WITH_INTEL="system"
WITH_GCC="no"

# Math Libraries (Intel MKL recommended)
MATH_MODE="mkl"
WITH_MKL="system"  # Use system MKL when MATH_MODE is mkl

# MPI Implementation (Intel MPI recommended, but OpenMPI/MPICH also supported)
MPI_MODE="intelmpi"
WITH_INTELMPI="system"
WITH_OPENMPI="no"
WITH_MPICH="no"

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="no"  # MKL provides ScaLAPACK
WITH_FFTW="no"       # MKL provides FFTW interface
WITH_LIBXC="install"
WITH_ELPA="install"

# Utility Libraries
WITH_CEREAL="install"
WITH_RAPIDJSON="install"

# Advanced Features (EXX calculations)
WITH_LIBRI="install"
WITH_LIBCOMM="install"

# Optional Features (MLALGO support)
WITH_LIBTORCH="no"
WITH_LIBNPY="no"
WITH_NEP="no"

# ELPA-GPU Support (uncomment and modify as needed)
# ENABLE_CUDA="yes"
# GPU_VERSION="75"  # Check your GPU compute capability
# export CUDA_PATH="/usr/local/cuda"

# ============================================================================
# Execution Mode Control
# ============================================================================
# Dry-run mode: Show what would be done without actually executing
DRY_RUN_MODE="no"   # Set to "yes" to enable dry-run mode

# Pack-run mode: Only check and install required packages
PACK_RUN_MODE="no"  # Set to "yes" to enable pack-run mode

# ============================================================================
# Intel Compiler and MPI Options
# ============================================================================
# Intel Compiler Version Selection
WITH_INTEL_CLASSIC="no"  # Set to "yes" to use classic Intel compilers (icc/icpc/ifort)
                         # Set to "no" to use new Intel compilers (icx/icpx/ifx)
                         # Classic compilers needed for AMD-CPU or GPU-version

WITH_IFX="yes"           # Set to "yes" to use new Fortran compiler ifx (default)
                        # Set to "no" to use traditional ifort
                        # Only applies when WITH_INTEL_CLASSIC="no"

# Intel MPI Version Selection  
INTELMPI_CLASSIC="no"   # Set to "yes" to use classic Intel MPI wrappers (mpiicc/mpiicpc/mpiifort)
                        # Set to "no" to use new Intel MPI wrappers (mpiicx/mpiicpx/mpiifx)
                        # Should match WITH_INTEL_CLASSIC setting for consistency
                        # Classic wrappers recommended for older Intel OneAPI versions


# ============================================================================
# Package Version Selection (main/alt versions)
# ============================================================================
# Choose between main (latest stable) and alt (alternative/legacy) versions
# Refer to scripts/package_versions.sh for specific version numbers

CMAKE_VERSION="main"        # main=3.31.7, alt=3.30.5
ELPA_VERSION="alt"         # main=2025.06.001, alt=2024.05.001 for intel oneapi<2024.2
LIBXC_VERSION="main"        # main=7.0.0, alt=6.2.2
# Optional Libraries
LIBTORCH_VERSION="main"     # main=2.1.2, alt=1.12.1 (use alt for older GLIBC)
# Note: main(2.1.2) version of LibTorch need glibc > 2.27
# Note: alt(1.12.1) version of LibTorch cannot support DeePMD-Torch for DPA

# Note: Intel toolchain uses MKL for math libraries (FFTW, ScaLAPACK)
# so OpenBLAS, FFTW, and ScaLAPACK version selections are not applicable

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain_new.sh \
  --with-intel="$WITH_INTEL" \
  --with-gcc="$WITH_GCC" \
  --math-mode="$MATH_MODE" \
  --mpi-mode="$MPI_MODE" \
  --with-mkl="$WITH_MKL" \
  --with-openmpi="$WITH_OPENMPI" \
  --with-mpich="$WITH_MPICH" \
  --with-intelmpi="$WITH_INTELMPI" \
  --with-cmake="$WITH_CMAKE" \
  --with-scalapack="$WITH_SCALAPACK" \
  --with-libxc="$WITH_LIBXC" \
  --with-fftw="$WITH_FFTW" \
  --with-elpa="$WITH_ELPA" \
  --with-cereal="$WITH_CEREAL" \
  --with-rapidjson="$WITH_RAPIDJSON" \
  --with-libtorch="$WITH_LIBTORCH" \
  --with-nep="$WITH_NEP" \
  --with-libnpy="$WITH_LIBNPY" \
  --with-libri="$WITH_LIBRI" \
  --with-libcomm="$WITH_LIBCOMM" \
  --with-intel-classic="$WITH_INTEL_CLASSIC" \
  --with-ifx="$WITH_IFX" \
  --with-intel-mpi-classic="$INTELMPI_CLASSIC" \
  --package-version cmake:"$CMAKE_VERSION" \
  --package-version elpa:"$ELPA_VERSION" \
  --package-version libxc:"$LIBXC_VERSION" \
  --package-version libtorch:"$LIBTORCH_VERSION" \
  ${DRY_RUN_MODE:+$([ "$DRY_RUN_MODE" = "yes" ] && echo "--dry-run")} \
  ${PACK_RUN_MODE:+$([ "$PACK_RUN_MODE" = "yes" ] && echo "--pack-run")} \
  ${ENABLE_CUDA:+--enable-cuda} \
  ${GPU_VERSION:+--gpu-ver="$GPU_VERSION"} \
  "$@" \
  | tee compile.log
