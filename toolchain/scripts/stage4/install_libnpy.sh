#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# libnpy is not need any complex setting
# Only problem is the installation from github.com
# Libnpy is under active development, you can check the latest version in github yourself

# Last Update in 2025-0504

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load LibNPY package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libnpy:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libnpy:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "libnpy" "$version_suffix"
dirname="libnpy-${libnpy_ver}"
filename="libnpy-${libnpy_ver}.tar.gz"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libnpy" ] && rm "${BUILDDIR}/setup_libnpy"

LIBNPY_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libnpy" in
    __INSTALL__)
        echo "==================== Installing LIBNPY ===================="
        pkg_install_dir="${INSTALLDIR}/$dirname"
        #pkg_install_dir="${HOME}/lib/libnpy/${libnpy_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        url="https://codeload.github.com/llohse/libnpy/tar.gz/v${libnpy_ver}"
        if verify_checksums "${install_lock_file}"; then
            echo "$dirname is already installed, skipping it."
        else
            if [ -f $filename ]; then
                echo "$filename is found"
            else
                # download from github.com and checksum
                echo "===> Notice: This version of Libnpy is downloaded in GitHub Release <==="
                download_pkg_from_url "${libnpy_sha256}" "${filename}" "${url}"
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
        LIBNPY_CFLAGS="-I'${pkg_install_dir}'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBNPY from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        # Find libnpy header file and derive package root directory
        libnpy_header_path="$(find_in_paths "npy.hpp" $INCLUDE_PATHS)"
        if [ "$libnpy_header_path" != "__FALSE__" ]; then
            # Derive pkg_install_dir from found header path
            # npy.hpp -> get include dir -> get parent dir
            libnpy_include_dir="$(dirname "$libnpy_header_path")"
            pkg_install_dir="$(dirname "$libnpy_include_dir")"
            echo "Found libnpy at: $pkg_install_dir"
            LIBNPY_CFLAGS="-I'${libnpy_include_dir}'"
        else
            report_error "Cannot find npy.hpp in system paths"
            exit 1
        fi
        ;;
    __DONTUSE__) ;;
    
    *)
        echo "==================== Linking LIBNPY to user paths ===================="
        pkg_install_dir="${with_libnpy}"
        check_dir "${pkg_install_dir}"
        LIBNPY_CFLAGS="-I'${pkg_install_dir}'"
        ;;
esac
if [ "$with_libnpy" != "__DONTUSE__" ]; then
    if [ "$with_libnpy" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_libnpy"
prepend_path CPATH "$pkg_install_dir/include"
export CPATH="${pkg_install_dir}/include":\${CPATH}
EOF
        cat "${BUILDDIR}/setup_libnpy" >> $SETUPFILE
    fi
    cat << EOF >> "${BUILDDIR}/setup_libnpy"
export LIBNPY_CFLAGS="${LIBNPY_CFLAGS}"
export LIBNPY_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_libnpy"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libnpy"
