#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# libcomm is not need any complex setting
# Only problem is the installation from github.com
# LibComm is under highly-active development, the git submodule installation is more recommended

# Last Update in 2025-0504
# other contributor: Peize Lin

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load LibComm package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libcomm:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libcomm:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "libcomm" "$version_suffix"
dirname="LibComm-${libcomm_ver}"
filename="LibComm-${libcomm_ver}.tar.gz"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libcomm" ] && rm "${BUILDDIR}/setup_libcomm"

libcomm_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libcomm" in
    __INSTALL__)
        echo "==================== Installing LIBCOMM ===================="
        pkg_install_dir="${INSTALLDIR}/$dirname"
        #pkg_install_dir="${HOME}/lib/libcomm/${libcomm_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        # url="https://github.com/abacusmodeling/LibComm/archive/refs/tags/v${libcomm_ver}.tar.gz"
        # url construction rules:
        # - Branch names (master, main, develop) without v prefix
        # - Version tags (e.g., 1.0.0) with v prefix
        if [[ "${libcomm_ver}" =~ ^(master|main|develop)$ ]]; then
            url="https://codeload.github.com/abacusmodeling/LibComm/tar.gz/${libcomm_ver}"
        else
            url="https://codeload.github.com/abacusmodeling/LibComm/tar.gz/v${libcomm_ver}"
        fi
        if verify_checksums "${install_lock_file}"; then
            echo "$dirname is already installed, skipping it."
        else
            if [ -f $filename ]; then
                echo "$filename is found"
            else
                # download from github.com and checksum
                echo "===> Notice: This version of LibComm is downloaded in GitHub master repository  <==="
                download_pkg_from_url "${libcomm_sha256}" "${filename}" "${url}"
            fi
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d $dirname ] && rm -rf $dirname
            tar -xzf $filename
            mkdir -p "${pkg_install_dir}"
            cp -r $dirname/* "${pkg_install_dir}/"
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
        LIBCOMM_CFLAGS="-I'${pkg_install_dir}'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBCOMM from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        # Find libcomm header file and derive package root directory
        libcomm_header_path="$(find_in_paths "Comm/Comm_Tools.h" $INCLUDE_PATHS)"
        if [ "$libcomm_header_path" != "__FALSE__" ]; then
            # Derive pkg_install_dir from found header path
            # Comm/Comm_Tools.h -> remove /Comm/Comm_Tools.h -> get include dir -> get parent dir
            libcomm_include_dir="$(dirname "$(dirname "$libcomm_header_path")")"
            pkg_install_dir="$(dirname "$libcomm_include_dir")"
            echo "Found libcomm at: $pkg_install_dir"
            LIBCOMM_CFLAGS="-I'${libcomm_include_dir}'"
        else
            report_error "Cannot find Comm/Comm_Tools.h in system paths"
            exit 1
        fi
        ;;
    __DONTUSE__) ;;
    
    *)
        echo "==================== Linking LIBCOMM to user paths ===================="
        pkg_install_dir="${with_libcomm}"
        check_dir "${pkg_install_dir}"
        LIBCOMM_CFLAGS="-I'${pkg_install_dir}'"
        ;;
esac
if [ "$with_libcomm" != "__DONTUSE__" ]; then
    if [ "$with_libcomm" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_libcomm"
prepend_path CPATH "$pkg_install_dir/include"
export CPATH="${pkg_install_dir}/include":\${CPATH}
EOF
        cat "${BUILDDIR}/setup_libcomm" >> $SETUPFILE
    fi
    cat << EOF >> "${BUILDDIR}/setup_libcomm"
export LIBCOMM_CFLAGS="${libcomm_CFLAGS}"
export LIBCOMM_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_libcomm"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libcomm"
