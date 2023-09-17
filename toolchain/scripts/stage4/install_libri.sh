#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# libri is not need any complex setting
# Only problem is the installation from github.com

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libri_ver="0.1.0"
libri_sha256="4721276e35b64e96f349df9899039159da0215ae0c1df94523c5c25fab3f7f89"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
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
    url="https://github.com/abacusmodeling/LibRI/archive/refs/tags/v${libri_ver}.tar.gz"
    filename="LibRI-${libri_ver}.tar.gz"
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
        checksum "$filename" "$libri_sha256"
        fi
        echo "Installing from scratch into ${pkg_install_dir}"
        [ -d $dirname ] && rm -rf $dirname
        tar -xzf $filename
        cp -r $dirname "${pkg_install_dir}/"
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
        ;;
    __SYSTEM__)
        echo "==================== CANNOT Finding LIBRI from system paths NOW ===================="
        recommend_offline_installation $filename $url
        # How to do it in libri? -- Zhaoqing in 2023/08/23
        # check_lib -lxcf03 "libxc"
        # check_lib -lxc "libxc"
        # add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
        # add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
        ;;
    __DONTUSE__) ;;
    
    *)
    echo "==================== Linking LIBRI to user paths ===================="
    check_dir "${pkg_install_dir}"
    LIBRI_CFLAGS="-I'${pkg_install_dir}'"
    ;;
esac
if [ "$with_libri" != "__DONTUSE__" ]; then
    if [ "$with_libri" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_libri"
prepend_path CPATH "$pkg_install_dir/include"
export CPATH="${pkg_install_dir}/include":${CPATH}
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
