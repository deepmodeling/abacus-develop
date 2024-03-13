#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# RAPIDJSON is not need any complex setting
# Only problem is the installation from github.com

# Last Update in 2024-0119

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

rapidjson_ver="1.1.0"
rapidjson_sha256="bf7ced29704a1e696fbccf2a2b4ea068e7774fa37f6d7dd4039d0787f8bed98e"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_rapidjson" ] && rm "${BUILDDIR}/setup_rapidjson"

RAPIDJSON_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_rapidjson" in
  __INSTALL__)
    echo "==================== Installing RAPIDJSON ===================="
    dirname="rapidjson-${rapidjson_ver}"
    pkg_install_dir="${INSTALLDIR}/$dirname"
    #pkg_install_dir="${HOME}/lib/rapidjson/${rapidjson_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    url="https://github.com/Tencent/rapidjson/archive/refs/tags/v${rapidjson_ver}.tar.gz"
    filename="rapidjson-${rapidjson_ver}.tar.gz"
    if verify_checksums "${install_lock_file}"; then
        echo "$dirname is already installed, skipping it."
    else
        if [ -f $filename ]; then
        echo "$filename is found"
        else
        # download from github.com and checksum
            echo "wget --quiet $url -O $filename"
            if ! wget --quiet $url -O $filename; then
            report_error "failed to download $url"
            recommend_offline_installation $filename $url
            fi
        # checksum
        checksum "$filename" "$rapidjson_sha256"
        fi
        echo "Installing from scratch into ${pkg_install_dir}"
        [ -d $dirname ] && rm -rf $dirname
        tar -xzf $filename
        mkdir -p "${pkg_install_dir}"
        cp -r $dirname/* "${pkg_install_dir}/"
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
        ;;
    __SYSTEM__)
        echo "==================== CANNOT Finding RAPIDJSON from system paths NOW ===================="
        recommend_offline_installation $filename $url
        # How to do it in rapidjson? -- Zhaoqing in 2023/08/23
        # check_lib -lxcf03 "libxc"
        # check_lib -lxc "libxc"
        # add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
        # add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
        ;;
    __DONTUSE__) ;;
    
    *)
    echo "==================== Linking RAPIDJSON to user paths ===================="
    check_dir "${pkg_install_dir}"
    RAPIDJSON_CFLAGS="-I'${pkg_install_dir}'"
    ;;
esac
if [ "$with_rapidjson" != "__DONTUSE__" ]; then
    if [ "$with_rapidjson" != "__SYSTEM__" ]; then
    # LibRI deps should find rapidjson include in CPATH
        cat << EOF > "${BUILDDIR}/setup_rapidjson"
prepend_path CPATH "$pkg_install_dir/include"
export CPATH="${pkg_install_dir}/include:"\${CPATH}
EOF
        cat "${BUILDDIR}/setup_rapidjson" >> $SETUPFILE
    fi
    cat << EOF >> "${BUILDDIR}/setup_rapidjson"
export RAPIDJSON_CFLAGS="${RAPIDJSON_CFLAGS}"
export RAPIDJSON_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_rapidjson"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "rapidjson"
