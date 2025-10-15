#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err
# Users can easily modify these parameters to customize the build

# Compiler Configuration
TOOLCHAIN_COMPILER="gnu"
WITH_GCC="install"
WITH_INTEL="no"

# Math Libraries
MATH_MODE="openblas"
WITH_OPENBLAS="install"

# MPI Implementation
WITH_OPENMPI="install"
WITH_4TH_OPENMPI="no"  # Set to "yes" for OpenMPI v4

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="install"
WITH_LIBXC="install"
WITH_FFTW="install"
WITH_ELPA="install"

# Utility Libraries
WITH_CEREAL="install"
WITH_RAPIDJSON="install"

# Optional Features (DeepKS support)
WITH_LIBTORCH="install"  # Set to "install" for DeepKS support
WITH_LIBNPY="install"    # Set to "install" for DeepKS support

# Advanced Features (EXX calculations)
WITH_LIBRI="install"     # Set to "install" for EXX calculations
WITH_LIBCOMM="install"   # Set to "install" for advanced communication

# GPU Support (uncomment and modify as needed)
# ENABLE_CUDA="yes"
# GPU_VERSION="75"  # Check your GPU compute capability
# export CUDA_PATH="/usr/local/cuda"

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain_new.sh \
  --with-gcc="$WITH_GCC" \
  --with-intel="$WITH_INTEL" \
  --math-mode="$MATH_MODE" \
  --with-openblas="$WITH_OPENBLAS" \
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
  --with-4th-openmpi="$WITH_4TH_OPENMPI" \
  ${ENABLE_CUDA:+--enable-cuda} \
  ${GPU_VERSION:+--gpu-ver="$GPU_VERSION"} \
  "$@" \
  | tee compile.log
