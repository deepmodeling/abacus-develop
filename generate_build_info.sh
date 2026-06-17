#!/bin/bash

# ==============================================================================
# generate_build_info.sh
#
# A script to generate build_info.h for Makefile-based builds.
# It attempts to detect system information, compiler flags, and library versions
# to provide a feature-rich build report similar to the CMake system.
# ==============================================================================

# set -e # Exit immediately if a command exits with a non-zero status.

if [ -z "$1" ]; then
    echo "Usage: $0 <output_file_path>"
    exit 1
fi
OUTPUT_FILE=$1

# --- Helper Functions ---

# Function to get version from a command's output
# Usage: get_version_from_command "<command>" "<regex>"
get_version_from_command() {
    local cmd="$1"
    local regex="$2"
    if command -v "$cmd" > /dev/null; then
        local output=$($cmd 2>&1)
        if [[ $output =~ $regex ]]; then
            echo "${BASH_REMATCH[1]}"
            return 0
        fi
    fi
    return 1
}

# Function to get version from a header file
# Usage: get_version_from_header "<header_path>" "<regex>"
get_version_from_header() {
    local header_path="$1"
    local regex="$2"
    if [ -f "$header_path" ]; then
        local line=$(grep -E "$regex" "$header_path" | head -n 1)
        if [[ $line =~ $regex ]]; then
            echo "${BASH_REMATCH[1]}"
            return 0
        fi
    fi
    return 1
}

# --- Main Detection Logic ---

# 1. Basic Info
PLATFORM_NAME="CPU"
BUILD_USER=$(whoami)
BUILD_HOST=$(hostname)
CXX_COMPILER_PATH="${CXX:-g++}"
CXX_COMPILER_VERSION=$($CXX_COMPILER_PATH --version 2>&1 | head -n 1)
CXX_FLAGS="${CXXFLAGS:-}"
LINKER_FLAGS="${LDFLAGS:-}"
CUDA_FLAGS="${CUDAFLAGS:-}"

# Detect Platform based on environment variables
if [ "${USE_ROCM}" == "ON" ]; then PLATFORM_NAME="CPU + AMD ROCm"; fi
if [ "${USE_CUDA}" == "ON" ]; then PLATFORM_NAME="CPU + NVIDIA CUDA"; fi
if [ "${USE_ELPA}" == "ON" ] && [ "${ENABLE_LCAO}" == "ON" ]; then PLATFORM_NAME="${PLATFORM_NAME} + ELPA"; fi

# 2. MPI
MPI_IMPLEMENTATION="no"
MPI_VERSION="no"
if [ "${ENABLE_MPI}" == "ON" ]; then
    MPI_IMPLEMENTATION="Unknown"
    MPI_COMPILER="${MPI_CXX_COMPILER:-mpicxx}"
    if command -v "$MPI_COMPILER" > /dev/null; then
        # Intel MPI (oneAPI)
        version=$(get_version_from_command "$MPI_COMPILER --version" "Intel\(R\) oneAPI DPC\+\+/C\+\+ Compiler ([0-9]+\.[0-9]+\.[0-9]+)")
        if [ -n "$version" ]; then MPI_IMPLEMENTATION="Intel MPI"; MPI_VERSION="$version"; fi
        # Intel MPI (classic)
        if [ -z "$version" ]; then
            version=$(get_version_from_command "$MPI_COMPILER --version" "icpc \(ICC\) ([0-9]+\.[0-9]+\.[0-9]+)")
            if [ -n "$version" ]; then MPI_IMPLEMENTATION="Intel MPI"; MPI_VERSION="$version"; fi
        fi
        # OpenMPI
        if [ -z "$version" ]; then
            version=$(get_version_from_command "$MPI_COMPILER --version" "Open MPI ([0-9]+\.[0-9]+\.[0-9]+)")
            if [ -n "$version" ]; then MPI_IMPLEMENTATION="OpenMPI"; MPI_VERSION="$version"; fi
        fi
        # MPICH
        if [ -z "$version" ]; then
            version=$(get_version_from_command "$MPI_COMPILER --version" "MPICH Version: ([0-9]+\.[0-9]+\.[0-9]+)")
            if [ -n "$version" ]; then MPI_IMPLEMENTATION="MPICH"; MPI_VERSION="$version"; fi
        fi
        # Fallback via mpirun
        if [ -z "$version" ] && command -v mpirun > /dev/null; then
            mr_out=$(mpirun --version 2>&1)
            if [[ $mr_out =~ Open\ MPI\ ([0-9]+\.[0-9]+\.[0-9]+) ]]; then
                MPI_IMPLEMENTATION="OpenMPI"; MPI_VERSION="${BASH_REMATCH[1]}"
            elif [[ $mr_out =~ MPICH\ Version:\ ([0-9]+\.[0-9]+\.[0-9]+) ]]; then
                MPI_IMPLEMENTATION="MPICH"; MPI_VERSION="${BASH_REMATCH[1]}"
            fi
        fi
    fi
    if [ "$MPI_VERSION" == "no" ]; then MPI_VERSION="yes (version unknown)"; else MPI_VERSION="yes (v$MPI_VERSION)"; fi
else
    MPI_IMPLEMENTATION="no"
    MPI_VERSION="no"
fi

# 3. OpenMP
OPENMP_VERSION="no"
if [[ "$CXX_FLAGS" == *"-fopenmp"* ]]; then
    OPENMP_VERSION="yes (version unknown)"
    # Try to get version from compiler
    if command -v "$CXX_COMPILER_PATH" > /dev/null; then
        version=$($CXX_COMPILER_PATH -dM -E - < /dev/null 2>&1 | grep '#define _OPENMP ' | awk '{print $3}')
        if [ -n "$version" ]; then
            # Convert YYYYMM to a more readable format if possible
            year=$((version / 10000))
            month=$((version % 10000))
            OPENMP_VERSION="yes (v${year}.${month})"
        fi
    fi
fi

# 4. MKL
MKL_SUPPORT="no"
if [ -n "$MKLROOT" ]; then
    MKL_SUPPORT="yes (version unknown)"
    version=$(get_version_from_header "$MKLROOT/include/mkl_version.h" "INTEL_MKL_VERSION ([0-9]+)")
    if [ -n "$version" ]; then
        major=$((version / 10000))
        minor=$(( (version % 10000) / 100 ))
        patch=$((version % 100))
        MKL_SUPPORT="yes (v${major}.${minor}.${patch})"
    fi
fi

# 5. Libxc
LIBXC_VERSION="no"
if [ "${ENABLE_LIBXC}" == "ON" ]; then
    LIBXC_VERSION="yes (version unknown)"
    # Try pkg-config first
    if command -v pkg-config > /dev/null; then
        version=$(pkg-config --modversion libxc 2>/dev/null)
        if [ -n "$version" ]; then LIBXC_VERSION="yes (v$version)"; fi
    fi
    # Fallback to header
    if [[ "$LIBXC_VERSION" == *"unknown"* ]]; then
        # Assuming standard include paths
        version=$(get_version_from_header "/usr/include/libxc/libxc.h" "LIBXC_VERSION_ \"([0-9]+\.[0-9]+\.[0-9]+)\"")
        if [ -n "$version" ]; then LIBXC_VERSION="yes (v$version)"; fi
    fi
fi

# 6. CUDA
CUDA_VERSION="no"
if [ "${USE_CUDA}" == "ON" ]; then
    CUDA_VERSION="yes (version unknown)"
    if [ -n "$CUDA_HOME" ]; then
        version=$(get_version_from_command "$CUDA_HOME/bin/nvcc --version" "release ([0-9]+\.[0-9]+)")
        if [ -n "$version" ]; then CUDA_VERSION="yes (v$version)"; fi
    elif command -v nvcc > /dev/null; then
        version=$(get_version_from_command "nvcc --version" "release ([0-9]+\.[0-9]+)")
        if [ -n "$version" ]; then CUDA_VERSION="yes (v$version)"; fi
    fi
fi

# 7. ROCm/HIP
ROCM_VERSION="no"
if [ "${USE_ROCM}" == "ON" ]; then
    ROCM_VERSION="yes (version unknown)"
    if [ -n "$ROCM_PATH" ]; then
        version=$(get_version_from_command "$ROCM_PATH/bin/hipcc --version" "HIP version: ([0-9]+\.[0-9]+\.[0-9]+)")
        if [ -n "$version" ]; then ROCM_VERSION="yes (v$version)"; fi
    elif command -v hipcc > /dev/null; then
        version=$(get_version_from_command "hipcc --version" "HIP version: ([0-9]+\.[0-9]+\.[0-9]+)")
        if [ -n "$version" ]; then ROCM_VERSION="yes (v$version)"; fi
    fi
fi

# 8. DeePMD-kit
DEEPMD_VERSION="no"
if [ -n "$DeePMD_DIR" ]; then
    DEEPMD_VERSION="yes (version unknown)"
    version=$(get_version_from_header "$DeePMD_DIR/include/deepmd/version.h" "global_install_prefix=\".*deepmd-kit-([0-9]+\.[0-9]+\.[0-9]+)\"")
    if [ -n "$version" ]; then DEEPMD_VERSION="yes (v$version)"; else DEEPMD_VERSION="yes (path: $DeePMD_DIR)"; fi
fi

# 9. ELPA
ELPA_VERSION="no"
if [ "${USE_ELPA}" == "ON" ] && [ -n "$ELPA_DIR" ]; then
    ELPA_VERSION="yes (version unknown)"
    if [ -f "$ELPA_DIR/bin/elpa2_print_version" ]; then
        version=$("$ELPA_DIR/bin/elpa2_print_version" 2>/dev/null)
        if [ -n "$version" ]; then ELPA_VERSION="yes ($version)"; fi
    else
        ELPA_VERSION="yes (path: $ELPA_DIR)"
    fi
fi

# 10. Cereal
CEREAL_VERSION="no"
if [ -n "$CEREAL_DIR" ]; then
    CEREAL_VERSION="yes (version unknown)"
    major=$(get_version_from_header "$CEREAL_DIR/include/cereal/version.hpp" "__CEREAL_VERSION_MAJOR[[:space:]]+([0-9]+)")
    minor=$(get_version_from_header "$CEREAL_DIR/include/cereal/version.hpp" "__CEREAL_VERSION_MINOR[[:space:]]+([0-9]+)")
    patch=$(get_version_from_header "$CEREAL_DIR/include/cereal/version.hpp" "__CEREAL_VERSION_PATCH[[:space:]]+([0-9]+)")
    if [ -n "$major" ] && [ -n "$minor" ] && [ -n "$patch" ]; then
        CEREAL_VERSION="yes (v${major}.${minor}.${patch})"
    else
        CEREAL_VERSION="yes (path: $CEREAL_DIR)"
    fi
fi

# 11. LibRI
LIBRI_VERSION="no"
if [ -n "$LIBRI_DIR" ]; then
    LIBRI_VERSION="yes (version unknown)"
    major=$(get_version_from_header "$LIBRI_DIR/include/RI/version.h" "__LIBRI_VERSION_MAJOR[[:space:]]+([0-9]+)")
    minor=$(get_version_from_header "$LIBRI_DIR/include/RI/version.h" "__LIBRI_VERSION_MINOR[[:space:]]+([0-9]+)")
    patch=$(get_version_from_header "$LIBRI_DIR/include/RI/version.h" "__LIBRI_VERSION_PATCH[[:space:]]+([0-9]+)")
    if [ -n "$major" ] && [ -n "$minor" ] && [ -n "$patch" ]; then
        LIBRI_VERSION="yes (v${major}.${minor}.${patch})"
    else
        LIBRI_VERSION="yes (path: $LIBRI_DIR)"
    fi
fi

# 12. LibComm
LIBCOMM_VERSION="no"
if [ -n "$LIBCOMM_DIR" ]; then
    LIBCOMM_VERSION="yes (path: $LIBCOMM_DIR)"
fi

# 13. FFTW (non-MKL)
FFTW3_VERSION="no"
if [ -z "$MKLROOT" ] && [ -n "$FFTW3_INCLUDE_DIR" ]; then
    FFTW3_VERSION="yes (version unknown)"
    hdr="$FFTW3_INCLUDE_DIR/fftw3.h"
    version=$(get_version_from_header "$hdr" "FFTW_VERSION\s+\"([^\"]+)\"")
    if [ -n "$version" ]; then FFTW3_VERSION="yes (v$version)"; fi
fi

# CUDA-aware MPI
CUDA_AWARE_MPI="no"
if [ "${USE_CUDA}" == "ON" ]; then
    if command -v ompi_info > /dev/null; then
        out=$(ompi_info --parsable --all 2>/dev/null)
        if echo "$out" | grep -q "mpi_built_with_cuda_support:value:true"; then
            CUDA_AWARE_MPI="yes"
        else
            CUDA_AWARE_MPI="no (or undetectable)"
        fi
    else
        CUDA_AWARE_MPI="no (or undetectable)"
    fi
fi
# --- Final File Generation ---

INPUT_FILE="source_io/build_info.h.in"

# Use sed to replace all placeholders with detected values
# Note the use of different delimiters (#) for paths to avoid conflicts with /
sed \
    -e "s#@ABACUS_PLATFORM_NAME@#$PLATFORM_NAME#g" \
    -e "s#@ABACUS_BUILD_TYPE@#${BUILD_TYPE:-Release}#g" \
    -e "s#@ABACUS_BUILD_USER@#$BUILD_USER#g" \
    -e "s#@ABACUS_BUILD_HOST@#$BUILD_HOST#g" \
    -e "s#@ABACUS_CXX_COMPILER_ID@#${CXX_COMPILER_ID:-Unknown}#g" \
    -e "s#@ABACUS_CXX_COMPILER_PATH@#$CXX_COMPILER_PATH#g" \
    -e "s#@ABACUS_CXX_COMPILER_VERSION@#$CXX_COMPILER_VERSION#g" \
    -e "s#@ABACUS_CXX_FLAGS@#$CXX_FLAGS#g" \
    -e "s#@ABACUS_LINKER_FLAGS@#$LINKER_FLAGS#g" \
    -e "s#@ABACUS_CUDA_FLAGS@#$CUDA_FLAGS#g" \
    -e "s#@ABACUS_MPI_IMPLEMENTATION@#$MPI_IMPLEMENTATION#g" \
    -e "s#@ABACUS_MPI_VERSION@#$MPI_VERSION#g" \
    -e "s#@ABACUS_CUDA_AWARE_MPI@#$CUDA_AWARE_MPI#g" \
    -e "s#@ABACUS_OPENMP_VERSION@#$OPENMP_VERSION#g" \
    -e "s#@ABACUS_MKL_SUPPORT@#$MKL_SUPPORT#g" \
    -e "s#@ABACUS_LIBXC_VERSION@#$LIBXC_VERSION#g" \
    -e "s#@ABACUS_FFTW_VERSION@#$FFTW3_VERSION#g" \
    -e "s#@ABACUS_CUDA_VERSION@#$CUDA_VERSION#g" \
    -e "s#@ABACUS_ROCM_VERSION@#$ROCM_VERSION#g" \
    -e "s#@ABACUS_DEEPMD_VERSION@#$DEEPMD_VERSION#g" \
    -e "s#@ABACUS_ELPA_VERSION@#$ELPA_VERSION#g" \
    -e "s#@ABACUS_CEREAL_VERSION@#$CEREAL_VERSION#g" \
    -e "s#@ABACUS_LIBRI_VERSION@#$LIBRI_VERSION#g" \
    -e "s#@ABACUS_LIBCOMM_VERSION@#$LIBCOMM_VERSION#g" \
    -e "s#@ABACUS_ASAN_STATUS@#${ENABLE_ASAN:-no}#g" \
    -e "s#@ABACUS_CMAKE_OPTIONS@#Not available with Makefile#g" \
    -e "s#@ABACUS_CMAKE_FIND_PACKAGES@#Not available with Makefile#g" \
    "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Generated $OUTPUT_FILE successfully."
