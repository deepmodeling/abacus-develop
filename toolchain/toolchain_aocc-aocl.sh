#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# Users can easily modify these parameters to customize the build
# Before running this script, ensure you have loaded your system packages

# Compiler Configuration
TOOLCHAIN_COMPILER="aocc-aocl"
WITH_AMD="system"
WITH_INTEL="no"
WITH_GCC="no"

# Math Libraries (AMD AOCL recommended)
MATH_MODE="aocl"
WITH_AOCL="system"

# MPI Implementation (OpenMPI recommended)
MPI_MODE="openmpi"
WITH_OPENMPI="install"
WITH_MPICH="no"
WITH_4TH_OPENMPI="no"  # Set to "yes" for OpenMPI v4, deprecated

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="system"  # AOCL provides ScaLAPACK
WITH_FFTW="system"       # AOCL provides FFTW
WITH_LIBXC="install"
WITH_ELPA="install"

# Utility Libraries
WITH_CEREAL="install"
WITH_RAPIDJSON="install"

# Advanced Features (EXX calculations)
WITH_LIBRI="install"
WITH_LIBCOMM="install"

# Optional Features (DeepKS support)
WITH_LIBTORCH="no"  # Set to "install" for DeepKS support
WITH_LIBNPY="no"    # Set to "install" for DeepKS support
WITH_NEP="no"

# AMD Compiler Options
WITH_FLANG="no"     # Flang is not recommended in current stage

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
ELPA_VERSION="main"         # main=2025.06.001, alt=2024.05.001
LIBXC_VERSION="main"        # main=7.0.0, alt=6.2.2
LIBTORCH_VERSION="main"     # main=2.1.2, alt=1.12.1 (use alt for older GLIBC)

# Note: AOCC-AOCL toolchain uses AOCL for math libraries (FFTW, ScaLAPACK)
# so OpenBLAS and ScaLAPACK version selections are not applicable

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain_new.sh \
  --with-amd="$WITH_AMD" \
  --with-gcc="$WITH_GCC" \
  --math-mode="$MATH_MODE" \
  --mpi-mode="$MPI_MODE" \
  --with-aocl="$WITH_AOCL" \
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
  --with-flang="$WITH_FLANG" \
  --with-4th-openmpi="$WITH_4TH_OPENMPI" \
  --package-version cmake:"$CMAKE_VERSION" \
  --package-version openmpi:"$OPENMPI_VERSION" \
  --package-version elpa:"$ELPA_VERSION" \
  --package-version libxc:"$LIBXC_VERSION" \
  --package-version libtorch:"$LIBTORCH_VERSION" \
  ${DRY_RUN_MODE:+$([ "$DRY_RUN_MODE" = "yes" ] && echo "--dry-run")} \
  ${PACK_RUN_MODE:+$([ "$PACK_RUN_MODE" = "yes" ] && echo "--pack-run")} \
  ${ENABLE_CUDA:+--enable-cuda} \
  ${GPU_VERSION:+--gpu-ver="$GPU_VERSION"} \
  "$@" \
  | tee compile.log
