#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# Users can easily modify these parameters to customize the build

# Compiler Configuration
TOOLCHAIN_COMPILER="intel"
WITH_INTEL="system"
WITH_GCC="no"

# Math Libraries (Intel MKL recommended)
MATH_MODE="mkl"

# MPI Implementation (Intel MPI recommended, but OpenMPI/MPICH also supported)
MPI_IMPLEMENTATION="intelmpi"
WITH_INTELMPI="system"

# Core Dependencies
WITH_CMAKE="install"
WITH_SCALAPACK="no"  # MKL provides ScaLAPACK
WITH_LIBXC="install"
WITH_FFTW="no"       # MKL provides FFTW interface
WITH_ELPA="install"

# Utility Libraries
WITH_CEREAL="install"
WITH_RAPIDJSON="install"

# Optional Features (DeepKS support)
WITH_LIBTORCH="no"  # Set to "install" for DeepKS support
WITH_LIBNPY="no"    # Set to "install" for DeepKS support

# Advanced Features (EXX calculations)
WITH_LIBRI="install"     # Recommended for EXX calculations
WITH_LIBCOMM="install"   # Recommended for advanced communication

# Intel Compiler Options
WITH_INTEL_CLASSIC="no"  # Set to "yes" for AMD-CPU or GPU-version

# GPU Support (uncomment and modify as needed)
# ENABLE_CUDA="yes"
# GPU_VERSION="75"  # Check your GPU compute capability
# export CUDA_PATH="/usr/local/cuda"

# ============================================================================
# System Requirements
# ============================================================================
# Before running this script, ensure you have loaded the required modules:
# module load mkl mpi compiler

# ============================================================================
# Execute Installation (DO NOT MODIFY BELOW THIS LINE)
# ============================================================================

# Call the main installation script with configured parameters
exec ./install_abacus_toolchain.sh \
  --with-intel="$WITH_INTEL" \
  --with-gcc="$WITH_GCC" \
  --math-mode="$MATH_MODE" \
  --with-mkl="$WITH_MKL" \
  --with-intelmpi="$WITH_INTELMPI" \
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
