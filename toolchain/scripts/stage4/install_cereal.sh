#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# CEREAL is not need any complex setting
# Only problem is the installation from github.com

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load Cereal package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "cereal:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "cereal:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "cereal" "$version_suffix"
dirname="cereal-${cereal_ver}"
filename="cereal-${cereal_ver}.tar.gz"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cereal" ] && rm "${BUILDDIR}/setup_cereal"

CEREAL_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_cereal" in
    __INSTALL__)
        echo "==================== Installing CEREAL ===================="
        pkg_install_dir="${INSTALLDIR}/$dirname"
        #pkg_install_dir="${HOME}/lib/cereal/${cereal_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        # url construction rules:
        # - Branch names (master, main, develop) without v prefix
        # - Version tags (e.g., 1.0.0) with v prefix
        if [[ "${cereal_ver}" =~ ^(master|main|develop)$ ]]; then
            url="https://codeload.github.com/USCiLab/cereal/tar.gz/${cereal_ver}"
        else
            url="https://codeload.github.com/USCiLab/cereal/tar.gz/v${cereal_ver}"
        fi
        if verify_checksums "${install_lock_file}"; then
            echo "$dirname is already installed, skipping it."
        else
            if [ -f $filename ]; then
                echo "$filename is found"
            else
                # download from github.com and checksum
                echo "===> Notice: This version of CEREAL is downloaded in GitHub master repository  <==="
                download_pkg_from_url "${cereal_sha256}" "${filename}" "${url}"
            fi
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d $dirname ] && rm -rf $dirname
            tar -xzf $filename
            #unzip -q $filename
            # apply patch files for libri installation in issue #6190, Kai Luo
            # echo ${SCRIPT_DIR}
            cd $dirname && pwd && patch -p1 < ${SCRIPT_DIR}/patches/6190.patch
            cd "${BUILDDIR}"
            # 
            mkdir -p "${pkg_install_dir}"
            cp -r $dirname/* "${pkg_install_dir}/"
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
        CEREAL_CFLAGS="-I'${pkg_install_dir}'"
        ;;
    __SYSTEM__)
        echo "==================== Finding CEREAL from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        # Find cereal header file and derive package root directory
        cereal_header_path="$(find_in_paths "cereal/cereal.hpp" $INCLUDE_PATHS)"
        if [ "$cereal_header_path" != "__FALSE__" ]; then
            # Derive pkg_install_dir from found header path
            # cereal/cereal.hpp -> remove /cereal/cereal.hpp -> get include dir -> get parent dir
            cereal_include_dir="$(dirname "$(dirname "$cereal_header_path")")"
            pkg_install_dir="$(dirname "$cereal_include_dir")"
            echo "Found cereal at: $pkg_install_dir"
            CEREAL_CFLAGS="-I'${cereal_include_dir}'"
        else
            report_error "Cannot find cereal/cereal.hpp in system paths"
            exit 1
        fi
        ;;
    __DONTUSE__) ;;
    
    *)
        echo "==================== Linking CEREAL to user paths ===================="
        pkg_install_dir="${with_cereal}"
        check_dir "${pkg_install_dir}"
        CEREAL_CFLAGS="-I'${pkg_install_dir}'"
        ;;
esac
if [ "$with_cereal" != "__DONTUSE__" ]; then
    if [ "$with_cereal" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_cereal"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}/include"
export CPATH="${pkg_install_dir}/include":\${CPATH}
export CMAKE_PREFIX_PATH="${pkg_install_dir}/include":\${CMAKE_PREFIX_PATH}
export CEREAL_ROOT="$pkg_install_dir"
EOF
    else
        cat << EOF > "${BUILDDIR}/setup_cereal"
export CEREAL_ROOT="$pkg_install_dir"
EOF
    fi
    cat "${BUILDDIR}/setup_cereal" >> $SETUPFILE
    cat << EOF >> "${BUILDDIR}/setup_cereal"
export CEREAL_CFLAGS="${CEREAL_CFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__CEREAL"
export CP_CFLAGS="\${CP_CFLAGS} ${CEREAL_CFLAGS}"
export CEREAL_VERSION="${cereal_ver}"
EOF
fi

load "${BUILDDIR}/setup_cereal"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cereal"
