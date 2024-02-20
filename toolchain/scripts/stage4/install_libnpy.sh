#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# libnpy is not need any complex setting
# Only problem is the installation from github.com
# Libnpy is under active development, you can check the latest version in github yourself

# Last Update in 2023-1124

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libnpy_ver="1.0.1"
libnpy_sha256="43452a4db1e8c1df606c64376ea1e32789124051d7640e7e4e8518ab4f0fba44"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libnpy" ] && rm "${BUILDDIR}/setup_libnpy"

LIBNPY_CFLAGS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libnpy" in
  __INSTALL__)
    echo "==================== Installing LIBNPY ===================="
    dirname="libnpy-${libnpy_ver}"
    pkg_install_dir="${INSTALLDIR}/$dirname"
    #pkg_install_dir="${HOME}/lib/libnpy/${libnpy_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    url="https://github.com/llohse/libnpy/archive/refs/tags/v${libnpy_ver}.tar.gz"
    filename="libnpy-${libnpy_ver}.tar.gz"
    if verify_checksums "${install_lock_file}"; then
        echo "$dirname is already installed, skipping it."
    else
        if [ -f $filename ]; then
        echo "$filename is found"
        else
        # download from github.com and checksum
            echo "===> Notice: This version Libnpy is downloaded in GitHub Release, which will always be out-of-date version <==="
            echo "wget --quiet $url -O $filename"
            if ! wget --quiet $url -O $filename; then
            report_error "failed to download $url"
            recommend_offline_installation $filename $url
            fi
        # checksum
        checksum "$filename" "$libnpy_sha256"
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
        echo "==================== CANNOT Finding LIBNPY from system paths NOW ===================="
        recommend_offline_installation $filename $url
        # How to do it in libnpy? -- Zhaoqing in 2023/08/23
        # check_lib -lxcf03 "libxc"
        # check_lib -lxc "libxc"
        # add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
        # add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
        ;;
    __DONTUSE__) ;;
    
    *)
    echo "==================== Linking LIBNPY to user paths ===================="
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
