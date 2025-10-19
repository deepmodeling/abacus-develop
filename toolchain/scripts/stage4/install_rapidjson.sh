#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# RapidJSON is not need any complex setting
# Only problem is the installation from github.com

# other contributor: Kai Luo, XingLiang Peng

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load RapidJSON package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "rapidjson:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "rapidjson:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "rapidjson" "$version_suffix"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_rapidjson" ] && rm "${BUILDDIR}/setup_rapidjson"

RAPIDJSON_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_rapidjson" in
    __INSTALL__)
        echo "==================== Installing RapidJSON ===================="
        dirname="rapidjson-${rapidjson_ver}"
        pkg_install_dir="${INSTALLDIR}/$dirname"
        #pkg_install_dir="${HOME}/lib/rapidjson/${rapidjson_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        # url construction rules:
        # - Branch names (master, main, develop) without v prefix
        # - Version tags (e.g., 1.0.0) with v prefix
        if [[ "${rapidjson_ver}" =~ ^(master|main|develop)$ ]]; then
            url="https://codeload.github.com/Tencent/rapidjson/tar.gz/${rapidjson_ver}"
        else
            url="https://codeload.github.com/Tencent/rapidjson/tar.gz/v${rapidjson_ver}"
        fi
        filename="rapidjson-${rapidjson_ver}.tar.gz"
        if verify_checksums "${install_lock_file}"; then
            echo "$dirname is already installed, skipping it."
        else
            if [ -f $filename ]; then
                echo "$filename is found"
            else
                # download from github.com and checksum
                echo "===> Notice: This version of RapidJSON is downloaded in GitHub master repository  <==="
                download_pkg_from_url "${rapidjson_sha256}" "${filename}" "${url}"
            fi
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d $dirname ] && rm -rf $dirname
            tar -xzf $filename
            #unzip -q $filename
            mkdir -p "${pkg_install_dir}"
            cp -r $dirname/* "${pkg_install_dir}/"
            # for rapidjson found in cmake
            cat << EOF > "${pkg_install_dir}/RapidJSONConfig.cmake"
get_filename_component(RAPIDJSON_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(RAPIDJSON_INCLUDE_DIRS "@INCLUDE_INSTALL_DIR@")
message(STATUS "RapidJSON found. Headers: ${RAPIDJSON_INCLUDE_DIRS}")
EOF
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
        RAPIDJSON_CFLAGS="-I'${pkg_install_dir}/include'"
        ;;
    __SYSTEM__)
        echo "==================== Finding RapidJSON from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        # Find rapidjson header file and derive package root directory
        rapidjson_header_path="$(find_in_paths "rapidjson/rapidjson.h" $INCLUDE_PATHS)"
        if [ "$rapidjson_header_path" != "__FALSE__" ]; then
            # Derive pkg_install_dir from found header path
            # rapidjson/rapidjson.h -> remove /rapidjson/rapidjson.h -> get include dir -> get parent dir
            rapidjson_include_dir="$(dirname "$(dirname "$rapidjson_header_path")")"
            pkg_install_dir="$(dirname "$rapidjson_include_dir")"
            echo "Found rapidjson at: $pkg_install_dir"
            RAPIDJSON_CFLAGS="-I'${rapidjson_include_dir}'"
        else
            report_error "Cannot find rapidjson/rapidjson.h in system paths"
            exit 1
        fi
        ;;
    __DONTUSE__) ;;

    *)
        echo "==================== Linking RapidJSON to user paths ===================="
        pkg_install_dir="${with_rapidjson}"
        check_dir "${pkg_install_dir}"
        RAPIDJSON_CFLAGS="-I'${pkg_install_dir}/include'"
        ;;
esac
if [ "$with_rapidjson" != "__DONTUSE__" ]; then
    if [ "$with_rapidjson" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_rapidjson"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export CPATH="$pkg_install_dir/include":\${CPATH}
export CMAKE_PREFIX_PATH="$pkg_install_dir":\${CMAKE_PREFIX_PATH}
export RAPIDJSON_ROOT="$pkg_install_dir"
EOF
    else
        cat << EOF > "${BUILDDIR}/setup_rapidjson"
export RAPIDJSON_ROOT="$pkg_install_dir"
EOF
    fi
    cat "${BUILDDIR}/setup_rapidjson" >> $SETUPFILE
    cat << EOF >> "${BUILDDIR}/setup_rapidjson"
export RAPIDJSON_CFLAGS="${RAPIDJSON_CFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__RAPIDJSON"
export CP_CFLAGS="\${CP_CFLAGS} ${RAPIDJSON_CFLAGS}"
export RAPIDJSON_VERSION="${rapidjson_ver}"
EOF
fi

load "${BUILDDIR}/setup_rapidjson"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "rapidjson"
