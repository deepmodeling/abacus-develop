#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# Users can easily modify these parameters to customize the build

# Compiler Configuration
TOOLCHAIN_COMPILER="gcc-aocl"
WITH_GCC="system"
WITH_AMD="no"
WITH_INTEL="no"

# Math Libraries (AMD AOCL recommended)
MATH_MODE="aocl"

# MPI Implementation (OpenMPI recommended)
WITH_OPENMPI="install"
WITH_4TH_OPENMPI="no"  # Set to "yes" for OpenMPI v4

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="system"  # AOCL provides ScaLAPACK
WITH_LIBXC="install"
WITH_FFTW="system"       # AOCL provides FFTW
WITH_ELPA="install"

# Utility Libraries
WITH_CEREAL="install"
WITH_RAPIDJSON="install"

# Optional Features (DeepKS support)
WITH_LIBTORCH="no"  # Set to "install" for DeepKS support
WITH_LIBNPY="no"    # Set to "install" for DeepKS support

# Advanced Features (EXX calculations)
WITH_LIBRI="no"     # Set to "install" for EXX calculations
WITH_LIBCOMM="no"   # Set to "install" for advanced communication

# GPU Support (uncomment and modify as needed)
# ENABLE_CUDA="yes"
# GPU_VERSION="75"  # Check your GPU compute capability
# export CUDA_PATH="/usr/local/cuda"

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain.sh \
  --with-gcc="$WITH_GCC" \
  --with-amd="$WITH_AMD" \
  --math-mode="$MATH_MODE" \
  --with-aocl="$WITH_AOCL" \
  --with-openmpi="$WITH_OPENMPI" \
  --with-cmake="$WITH_CMAKE" \
  --with-scalapack="$WITH_SCALAPACK" \
  --with-libxc="$WITH_LIBXC" \
  --with-fftw="$WITH_FFTW" \
  --with-elpa="$WITH_ELPA" \
  --with-cereal="$WITH_CEREAL" \
  --with-rapidjson="$WITH_RAPIDJSON" \
  --with-libtorch="$WITH_LIBTORCH" \
  --with-libnpy="$WITH_LIBNPY" \
  --with-libri="$WITH_LIBRI" \
  --with-libcomm="$WITH_LIBCOMM" \
  ${ENABLE_CUDA:+--enable-cuda} \
  ${GPU_VERSION:+--gpu-ver="$GPU_VERSION"} \
  "$@" \
  | tee compile.log