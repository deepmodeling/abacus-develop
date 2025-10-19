#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err
# Users can easily modify these parameters to customize the build

# Before running this script, ensure you have loaded your system packages

# Compiler Configuration
TOOLCHAIN_COMPILER="gnu"
WITH_GCC="system"
WITH_INTEL="no"
WITH_AMD="no"

# Math Libraries
MATH_MODE="openblas"
WITH_OPENBLAS="install"

# MPI Implementation (OpenMPI recommended)
MPI_MODE="openmpi"
WITH_OPENMPI="install"
WITH_4TH_OPENMPI="no"  # Set to "yes" for OpenMPI v4, deprecated
WITH_MPICH="no"

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="install"
WITH_LIBXC="install"
WITH_FFTW="install"
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
# Package Version Selection (main/alt versions)
# ============================================================================
# Choose between main (latest stable) and alt (alternative/legacy) versions
# Refer to scripts/package_versions.sh for specific version numbers

CMAKE_VERSION="main"        # main=3.31.7, alt=3.30.5
OPENMPI_VERSION="main"      # main=5.0.8, alt=4.1.6
OPENBLAS_VERSION="main"     # main=0.3.30, alt=0.3.27
ELPA_VERSION="main"         # main=2025.06.001, alt=2024.05.001
LIBXC_VERSION="main"        # main=7.0.0, alt=6.2.2
SCALAPACK_VERSION="main"    # main=2.2.2, alt=2.2.1
# Optional Libraries
LIBTORCH_VERSION="main"     # main=2.1.2, alt=1.12.1 (use alt for older GLIBC)

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain_new.sh \
  --with-gcc="$WITH_GCC" \
  --with-intel="$WITH_INTEL" \
  --with-amd="$WITH_AMD" \
  --math-mode="$MATH_MODE" \
  --mpi-mode="$MPI_MODE" \
  --with-openblas="$WITH_OPENBLAS" \
  --with-openmpi="$WITH_OPENMPI" \
  --with-mpich="$WITH_MPICH" \
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
  --with-4th-openmpi="$WITH_4TH_OPENMPI" \
  --package-version cmake:"$CMAKE_VERSION" \
  --package-version openmpi:"$OPENMPI_VERSION" \
  --package-version openblas:"$OPENBLAS_VERSION" \
  --package-version elpa:"$ELPA_VERSION" \
  --package-version libxc:"$LIBXC_VERSION" \
  --package-version scalapack:"$SCALAPACK_VERSION" \
  --package-version libtorch:"$LIBTORCH_VERSION" \
  ${DRY_RUN_MODE:+$([ "$DRY_RUN_MODE" = "yes" ] && echo "--dry-run")} \
  ${PACK_RUN_MODE:+$([ "$PACK_RUN_MODE" = "yes" ] && echo "--pack-run")} \
  ${ENABLE_CUDA:+--enable-cuda} \
  ${GPU_VERSION:+--gpu-ver="$GPU_VERSION"} \
  "$@" \
  | tee compile.log
