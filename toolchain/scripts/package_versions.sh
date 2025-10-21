#!/bin/bash
# ABACUS Toolchain - Centralized Package Version Management
# 
# This file contains all package versions, checksums, and URLs in one place.
#
# Usage: source this file in install scripts to get package version information
# Example: source "${SCRIPT_DIR}/package_versions.sh" && load_package_vars "openmpi"

# =============================================================================
# STAGE 0: Build Tools and Compilers
# =============================================================================

# GCC (supports dual versions) - Special case: main=13.2.0, alt=11.4.0
gcc_main_ver="13.2.0"
gcc_main_sha256="8cb4be3796651976f94b9356fa08d833524f62420d6292c5033a9a26af315078"
gcc_alt_ver="11.4.0"
gcc_alt_sha256="af828619dd1970734dda3cfb792ea3f2cba61b5a00170ba8bce4910749d73c07"

# CMake (supports dual versions) - main=3.31.7, alt=3.30.5
cmake_main_ver="3.31.7"
cmake_main_sha256_x86_64="b7a5c909cdafc36042c8c9bd5765e92ff1f2528cf01720aa6dc4df294ec7e1a0"
cmake_main_sha256_aarch64="ce8e32b2c1c497dd7f619124c043ac5c28a88677e390c58748dd62fe460c62a2"
cmake_main_sha256_macos="1cb11aa2edae8551bb0f22807c6f5246bd0eb60ae9fa1474781eb4095d299aca"
cmake_alt_ver="3.30.5"
cmake_alt_sha256_x86_64="83de8839f3fb0d9caf982a0435da3fa8c4fbe2c817dfec99def310dc7e6a8404"
cmake_alt_sha256_aarch64="93c3b8920379585dece1314f113c6c9008eaedfe56023c78d856fc86dad5b8e2"
cmake_alt_sha256_macos="3d603e507c7579b13518ef752b4ffcf3ed479fba80ee171d7d85da8153e869d0"

# =============================================================================
# STAGE 1: MPI Implementations
# =============================================================================

# OpenMPI (supports dual versions) - main=5.0.8, alt=4.1.6
openmpi_main_ver="5.0.8"
openmpi_main_sha256="53131e1a57e7270f645707f8b0b65ba56048f5b5ac3f68faabed3eb0d710e449"
openmpi_alt_ver="4.1.6"
openmpi_alt_sha256="f740994485516deb63b5311af122c265179f5328a0d857a567b85db00b11e415"

# MPICH (supports dual versions) - main=4.3.1, alt=4.1.0
mpich_main_ver="4.3.1"
mpich_main_sha256="acc11cb2bdc69678dc8bba747c24a28233c58596f81f03785bf2b7bb7a0ef7dc"
mpich_alt_ver="4.1.0"
mpich_alt_sha256="8b1ec63bc44c7caa2afbb457bc5b3cd4a70dbe46baba700123d67c48dc5ab6a0"

# =============================================================================
# STAGE 2: Math Libraries
# =============================================================================

# OpenBLAS (supports dual versions) - main=0.3.30, alt=0.3.27
openblas_main_ver="0.3.30"
openblas_main_sha256="27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d"
openblas_alt_ver="0.3.27"
openblas_alt_sha256="aa2d68b1564fe2b13bc292672608e9cdeeeb6dc34995512e65c3b10f4599e897"

# =============================================================================
# STAGE 3: Scientific Computing Libraries
# =============================================================================

# ELPA (supports dual versions) - main=2025.06.001, alt=2024.05.001
elpa_main_ver="2025.06.001"
elpa_main_sha256="feeb1fea1ab4a8670b8d3240765ef0ada828062ef7ec9b735eecba2848515c94"
elpa_alt_ver="2024.05.001"
elpa_alt_sha256="9caf41a3e600e2f6f4ce1931bd54185179dade9c171556d0c9b41bbc6940f2f6"

# FFTW (supports dual versions) - Special case: both main and alt are 3.3.10
fftw_main_ver="3.3.10"
fftw_main_sha256="56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467"
fftw_alt_ver="3.3.10"
fftw_alt_sha256="56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467"

# LibXC (supports dual versions) - main=7.0.0, alt=6.2.2
libxc_main_ver="7.0.0"
libxc_main_sha256="e9ae69f8966d8de6b7585abd9fab588794ada1fab8f689337959a35abbf9527d"
libxc_alt_ver="6.2.2"
libxc_alt_sha256="f72ed08af7b9dff5f57482c5f97bff22c7dc49da9564bc93871997cbda6dacf3"

# ScaLAPACK (supports dual versions) - main=2.2.2, alt=2.2.1
scalapack_main_ver="2.2.2"
scalapack_main_sha256="a2f0c9180a210bf7ffe126c9cb81099cf337da1a7120ddb4cbe4894eb7b7d022"
scalapack_alt_ver="2.2.1"
scalapack_alt_sha256="4aede775fdb28fa44b331875730bcd5bab130caaec225fadeccf424c8fcb55aa"

# =============================================================================
# STAGE 4: Advanced Feature Libraries
# =============================================================================

# LibTorch (supports dual versions) - main=2.1.2, alt=1.12.1
libtorch_main_ver="2.1.2"
libtorch_main_sha256="904b764df6106a8a35bef64c4b55b8c1590ad9d071eb276e680cf42abafe79e9"
libtorch_alt_ver="1.12.1"
libtorch_alt_sha256="82c7be80860f2aa7963f8700004a40af8205e1d721298f2e09b700e766a9d283"
# user can manually download higher version of libtorch by:
# wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-{libtorch_ver}%2Bcpu.zip
# 2.1.2 recommended for lower GLIBC support (lower than 3.4.26)

# LibNPY (supports dual versions) - Special case: both main and alt are 1.0.1
libnpy_main_ver="1.0.1"
libnpy_main_sha256="43452a4db1e8c1df606c64376ea1e32789124051d7640e7e4e8518ab4f0fba44"
libnpy_alt_ver="1.0.1"
libnpy_alt_sha256="43452a4db1e8c1df606c64376ea1e32789124051d7640e7e4e8518ab4f0fba44"

# Master branch packages (no fixed versions)
cereal_ver="master"
cereal_sha256="--no-checksum"

libcomm_ver="master"
libcomm_sha256="--no-checksum"

libri_ver="master"
libri_sha256="--no-checksum"

rapidjson_ver="master"
rapidjson_sha256="--no-checksum"

# NEP (Neural Evolution Potential) - CPU version
nep_ver="main"
nep_sha256="--no-checksum"

# =============================================================================
# Package Variable Loading Function
# =============================================================================

load_package_vars() {
    local package_name="$1"
    local version_suffix="$2"  # Optional version suffix for multi-version packages
    
    case "${package_name}" in
        "gcc")
            if [ "${version_suffix}" = "alt" ]; then
                gcc_ver="${gcc_alt_ver}"
                gcc_sha256="${gcc_alt_sha256}"
            else
                gcc_ver="${gcc_main_ver}"
                gcc_sha256="${gcc_main_sha256}"
            fi
            ;;
        "cmake")
            # Determine architecture for SHA256 selection
            local arch_suffix=""
            if [ "${OPENBLAS_ARCH}" = "arm64" ]; then
                if [ "$(uname -s)" = "Darwin" ]; then
                    arch_suffix="_macos"
                else
                    arch_suffix="_aarch64"
                fi
            else
                arch_suffix="_x86_64"
            fi
            
            if [ "${version_suffix}" = "alt" ]; then
                cmake_ver="${cmake_alt_ver}"
                eval "cmake_sha256=\${cmake_alt_sha256${arch_suffix}}"
            else
                cmake_ver="${cmake_main_ver}"
                eval "cmake_sha256=\${cmake_main_sha256${arch_suffix}}"
            fi
            ;;
        "openmpi")
            if [ "${OPENMPI_4TH}" = "yes" ]; then
                echo "WARNING: OPENMPI_4TH=yes is deprecated. Please use 'alt' parameter instead." >&2
                openmpi_ver="${openmpi_alt_ver}"
                openmpi_sha256="${openmpi_alt_sha256}"
            elif [ "${version_suffix}" = "alt" ]; then
                openmpi_ver="${openmpi_alt_ver}"
                openmpi_sha256="${openmpi_alt_sha256}"
            else
                openmpi_ver="${openmpi_main_ver}"
                openmpi_sha256="${openmpi_main_sha256}"
            fi
            ;;
        "mpich")
            if [ "${version_suffix}" = "alt" ]; then
                mpich_ver="${mpich_alt_ver}"
                mpich_sha256="${mpich_alt_sha256}"
            else
                mpich_ver="${mpich_main_ver}"
                mpich_sha256="${mpich_main_sha256}"
            fi
            ;;
        "openblas")
            if [ "${version_suffix}" = "alt" ]; then
                openblas_ver="${openblas_alt_ver}"
                openblas_sha256="${openblas_alt_sha256}"
            else
                openblas_ver="${openblas_main_ver}"
                openblas_sha256="${openblas_main_sha256}"
            fi
            ;;
        "elpa")
            if [ "${version_suffix}" = "alt" ]; then
                elpa_ver="${elpa_alt_ver}"
                elpa_sha256="${elpa_alt_sha256}"
            else
                elpa_ver="${elpa_main_ver}"
                elpa_sha256="${elpa_main_sha256}"
            fi
            ;;
        "fftw")
            if [ "${version_suffix}" = "alt" ]; then
                fftw_ver="${fftw_alt_ver}"
                fftw_sha256="${fftw_alt_sha256}"
            else
                fftw_ver="${fftw_main_ver}"
                fftw_sha256="${fftw_main_sha256}"
            fi
            ;;
        "libxc")
            if [ "${version_suffix}" = "alt" ]; then
                libxc_ver="${libxc_alt_ver}"
                libxc_sha256="${libxc_alt_sha256}"
            else
                libxc_ver="${libxc_main_ver}"
                libxc_sha256="${libxc_main_sha256}"
            fi
            ;;
        "scalapack")
            if [ "${version_suffix}" = "alt" ]; then
                scalapack_ver="${scalapack_alt_ver}"
                scalapack_sha256="${scalapack_alt_sha256}"
            else
                scalapack_ver="${scalapack_main_ver}"
                scalapack_sha256="${scalapack_main_sha256}"
            fi
            ;;
        "libtorch")
            if [ "${version_suffix}" = "alt" ]; then
                libtorch_ver="${libtorch_alt_ver}"
                libtorch_sha256="${libtorch_alt_sha256}"
            else
                libtorch_ver="${libtorch_main_ver}"
                libtorch_sha256="${libtorch_main_sha256}"
            fi
            ;;
        "libnpy")
            if [ "${version_suffix}" = "alt" ]; then
                libnpy_ver="${libnpy_alt_ver}"
                libnpy_sha256="${libnpy_alt_sha256}"
            else
                libnpy_ver="${libnpy_main_ver}"
                libnpy_sha256="${libnpy_main_sha256}"
            fi
            ;;
        "cereal")
            cereal_ver="${cereal_ver}"
            cereal_sha256="${cereal_sha256}"
            ;;
        "libcomm")
            libcomm_ver="${libcomm_ver}"
            libcomm_sha256="${libcomm_sha256}"
            ;;
        "libri")
            libri_ver="${libri_ver}"
            libri_sha256="${libri_sha256}"
            ;;
        "rapidjson")
            rapidjson_ver="${rapidjson_ver}"
            rapidjson_sha256="${rapidjson_sha256}"
            ;;
        "nep")
            nep_ver="${nep_ver}"
            nep_sha256="${nep_sha256}"
            ;;
        *)
            echo "Error: Unknown package '${package_name}'"
            return 1
            ;;
    esac
}