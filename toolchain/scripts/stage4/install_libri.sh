#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# libri is not need any complex setting
# Only problem is the installation from github.com
# LibRI is under highly-active development, the git submodule installation is more recommended

# Last Update in 2025-0504
# other contributor: Peize Lin

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load LibRI package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libri:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libri:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "libri" "$version_suffix"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libri" ] && rm "${BUILDDIR}/setup_libri"

libri_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libri" in
    __INSTALL__)
        echo "==================== Installing LIBRI ===================="
        dirname="LibRI-${libri_ver}"
        pkg_install_dir="${INSTALLDIR}/$dirname"
        #pkg_install_dir="${HOME}/lib/libri/${libri_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        # url construction rules:
        # - Branch names (master, main, develop) without v prefix
        # - Version tags (e.g., 1.0.0) with v prefix
        if [[ "${libri_ver}" =~ ^(master|main|develop)$ ]]; then
            url="https://codeload.github.com/abacusmodeling/LibRI/tar.gz/${libri_ver}"
        else
            url="https://codeload.github.com/abacusmodeling/LibRI/tar.gz/v${libri_ver}"
        fi
        filename="LibRI-${libri_ver}.tar.gz"
        if verify_checksums "${install_lock_file}"; then
            echo "$dirname is already installed, skipping it."
        else
            if [ -f $filename ]; then
                echo "$filename is found"
            else
                # download from github.com and checksum
                echo "===> Notice: This version of LibRI is downloaded in GitHub Release <==="
                download_pkg_from_url "${libri_sha256}" "${filename}" "${url}"
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
        LIBRI_CFLAGS="-I'${pkg_install_dir}'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBRI from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        # Find libri header file and derive package root directory
        libri_header_path="$(find_in_paths "RI/version.h" $INCLUDE_PATHS)"
        if [ "$libri_header_path" != "__FALSE__" ]; then
            # Derive pkg_install_dir from found header path
            # RI/version.h -> remove /RI/version.h -> get include dir -> get parent dir
            libri_include_dir="$(dirname "$(dirname "$libri_header_path")")"
            pkg_install_dir="$(dirname "$libri_include_dir")"
            echo "Found libri at: $pkg_install_dir"
            LIBRI_CFLAGS="-I'${libri_include_dir}'"
        else
            report_error "Cannot find RI/version.h in system paths"
            exit 1
        fi
        ;;
    __DONTUSE__) ;;

    *)
        echo "==================== Linking LIBRI to user paths ===================="
        pkg_install_dir="${with_libri}"
        check_dir "${pkg_install_dir}"
        LIBRI_CFLAGS="-I'${pkg_install_dir}'"
        ;;
esac
if [ "$with_libri" != "__DONTUSE__" ]; then
    if [ "$with_libri" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_libri"
prepend_path CPATH "$pkg_install_dir/include"
export CPATH="${pkg_install_dir}/include":\${CPATH}
EOF
        cat "${BUILDDIR}/setup_libri" >> $SETUPFILE
    fi
    cat << EOF >> "${BUILDDIR}/setup_libri"
export LIBRI_CFLAGS="${libri_CFLAGS}"
export LIBRI_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_libri"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libri"
