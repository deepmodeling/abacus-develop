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
cmake_main_sha256_x86_64="14e15d0b445dbeac686acc13fe13b3135e8307f69ccf4c5c91403996ce5aa2d4"
cmake_main_sha256_aarch64="e5b2dc2aefdca10afe09c8fa4ee2bbb4e732665943a94322f99c118781910c3c"
cmake_main_sha256_macos="1cb11aa2edae8551bb0f22807c6f5246bd0eb60ae9fa1474781eb4095d299aca"
cmake_alt_ver="3.30.5"
cmake_alt_sha256_x86_64="f747d9b23e1a252a8beafb4ed2bc2ddf78cff7f04a8e4de19f4ff88e9b51dc9d"
cmake_alt_sha256_aarch64="da7dead2c92c1747b40d506d7f7d68590f5bab175316d2e7af73e48a2e417e48"
cmake_alt_sha256_macos="3d603e507c7579b13518ef752b4ffcf3ed479fba80ee171d7d85da8153e869d0"

# =============================================================================
# STAGE 1: MPI Implementations
# =============================================================================

# OpenMPI (supports dual versions) - main=5.0.10, alt=4.1.8
openmpi_main_ver="5.0.10"
openmpi_main_sha256="0acecc4fc218e5debdbcb8a41d182c6b0f1d29393015ed763b2a91d5d7374cc6"
openmpi_alt_ver="4.1.8"
openmpi_alt_sha256="466f68e3132a1dc02710cc2011fafced8336d98359fa2dae4dddcfd5719f12a9"

# MPICH (supports dual versions) - main=5.0.1, alt=4.3.2
mpich_main_ver="5.0.1"
mpich_main_sha256="8c1832a13ddacf071685069f5fadfd1f2877a29e1a628652892c65211b1f3327"
mpich_alt_ver="4.3.2"
mpich_alt_sha256="47d774587a7156a53752218c811c852e70ac44db9c502dc3f399b4cb817e3818"

# =============================================================================
# STAGE 2: Math Libraries
# =============================================================================

# OpenBLAS (supports dual versions) - main=0.3.33, alt=0.3.30
openblas_main_ver="0.3.33"
openblas_main_sha256="6761af1d9f5d353ab4f0b7497be2643313b36c8f31caec0144bfef198e71e6ab"
openblas_alt_ver="0.3.30"
openblas_alt_sha256="27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d"

# =============================================================================
# STAGE 3: Scientific Computing Libraries
# =============================================================================

# ELPA (supports dual versions) - main=2026.02.001, alt=2024.05.001
elpa_main_ver="2026.02.001"
elpa_main_sha256="a379f27f4dbd27b2ee45017afec656d064301e97150c874649bdfd64957b75ed"
elpa_alt_ver="2024.05.001"
elpa_alt_sha256="9caf41a3e600e2f6f4ce1931bd54185179dade9c171556d0c9b41bbc6940f2f6"

# FFTW (supports dual versions) - main=3.3.11, alt=3.3.10
fftw_main_ver="3.3.11"
fftw_main_sha256="5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1"
fftw_alt_ver="3.3.10"
fftw_alt_sha256="56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467"

# LibXC (supports dual versions) - main=7.0.0, alt=6.2.2
libxc_main_ver="7.0.0"
libxc_main_sha256="e9ae69f8966d8de6b7585abd9fab588794ada1fab8f689337959a35abbf9527d"
libxc_alt_ver="6.2.2"
libxc_alt_sha256="f72ed08af7b9dff5f57482c5f97bff22c7dc49da9564bc93871997cbda6dacf3"

# ScaLAPACK (supports dual versions) - main=2.2.2, alt=2.2.1
scalapack_main_ver="2.2.3"
scalapack_main_sha256="5d93701eca663925e98010dd8d0f45fd79b2191d74e5afa5711d587370a8b9dd"
scalapack_alt_ver="2.2.2"
scalapack_alt_sha256="a2f0c9180a210bf7ffe126c9cb81099cf337da1a7120ddb4cbe4894eb7b7d022"

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

# Branch packages cut with specific commits
cereal_ver="22a1b36"
cereal_sha256="a8171736e6b553dd6cd37919c13433b01f499d24d45af502975a9439728803e0"

libcomm_ver="965bf90"
libcomm_sha256="d7b991465d98d7b715b484d86880bf3525b9bf0cc62c3e5d38b0f6d140f6b9d4"

libri_ver="e6d78e0"
libri_sha256="619b49a14047d7a167515d1f1d0fa2d82fbebd63b8cbd3181e07df6ed993a22c"

rapidjson_ver="24b5e7a"
rapidjson_sha256="dcb57b11036cb8fc6b2a57a6aded68d52e9cfe543811bf4fa8941087f84e72d0"

# NEP (Neural Evolution Potential) - CPU version
nep_ver="629ec5d"
nep_sha256="4d4d3c64211a2a39e5d5c795b77befbba987cc809786e0cd6abfa46d0f3bf8cb"

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
