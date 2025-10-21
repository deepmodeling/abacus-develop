#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

# Load centralized version management
source "${SCRIPT_DIR}"/package_versions.sh

# Load OpenBLAS package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "openblas:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "openblas:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "openblas" "$version_suffix"
openblas_pkg="OpenBLAS-${openblas_ver}.tar.gz"

find_openblas_dir() {
    local __dir=''
    for __dir in *OpenBLAS*; do
        if [ -d "$__dir" ]; then
            echo "$__dir"
            return
        fi
    done
    echo ''
}

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

echo "==================== Getting proc arch info using OpenBLAS tools ===================="
if [ -f ${openblas_pkg} ]; then
    echo "${openblas_pkg} is found"
else
# using codeload.github
    url="https://codeload.github.com/OpenMathLib/OpenBLAS/tar.gz/v${openblas_ver}"
    download_pkg_from_url "${openblas_sha256}" "${openblas_pkg}" "${url}"
fi
# if toolchain run in pack-run mode, do exit
if [ "${PACK_RUN}" = "__TRUE__" ]; then
    echo "--pack-run mode specified, skip arch detection"
    exit 0
fi
tar -xzf ${openblas_pkg}
openblas_dir="$(find_openblas_dir)"
openblas_conf="${openblas_dir}/Makefile.conf"
# try find Makefile.config, if not then generate one with make lapack_prebuild
if ! [ -f "$openblas_conf" ]; then
    cd "$openblas_dir"
    make lapack_prebuild
    cd ..
fi
OPENBLAS_LIBCORE="$(grep 'LIBCORE=' $openblas_conf | cut -f2 -d=)"
OPENBLAS_ARCH="$(grep 'ARCH=' $openblas_conf | cut -f2 -d=)"
echo "OpenBLAS detected LIBCORE = $OPENBLAS_LIBCORE"
echo "OpenBLAS detected ARCH    = $OPENBLAS_ARCH"
# output setup file
cat << EOF > "${BUILDDIR}/openblas_arch"
export OPENBLAS_LIBCORE="${OPENBLAS_LIBCORE}"
export OPENBLAS_ARCH="${OPENBLAS_ARCH}"
EOF
