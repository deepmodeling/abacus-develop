#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-0309

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# libxc_ver="6.2.2"
# libxc_sha256="a0f6f1bba7ba5c0c85b2bfe65aca1591025f509a7f11471b4cd651a79491b045"
libxc_ver="7.0.0"
libxc_sha256="e9ae69f8966d8de6b7585abd9fab588794ada1fab8f689337959a35abbf9527d"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxc" ] && rm "${BUILDDIR}/setup_libxc"

LIBXC_CFLAGS=""
LIBXC_LDFLAGS=""
LIBXC_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxc" in
  __INSTALL__)
    echo "==================== Installing LIBXC ===================="
    pkg_install_dir="${INSTALLDIR}/libxc-${libxc_ver}"
    #pkg_install_dir="${HOME}/lib/libxc/${libxc_ver}-gcc8"
    install_lock_file="$pkg_install_dir/install_successful"
    libxc_pkg="libxc-${libxc_ver}.tar.bz2"
    if verify_checksums "${install_lock_file}"; then
      echo "libxc-${libxc_ver} is already installed, skipping it."
    else
      if [ -f ${libxc_pkg} ]; then
        echo "${libxc_pkg} is found"
      else
        #download_pkg_from_ABACUS_org "${libxc_sha256}" "${libxc_pkg}"
        libxc_url="https://gitlab.com/libxc/libxc/-/archive/${libxc_ver}/${libxc_pkg}"
        download_pkg_from_url  "${libxc_sha256}" "${libxc_pkg}" "${libxc_url}"
      fi
    if [ "${PACK_RUN}" = "__TRUE__" ]; then
      echo "--pack-run mode specified, skip installation"
    else
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxc-${libxc_ver} ] && rm -rf libxc-${libxc_ver}
      tar -xjf ${libxc_pkg}
      cd libxc-${libxc_ver}
      # using cmake method to install libxc is neccessary for abacus
      mkdir build 
      cd build
      cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        -DBUILD_SHARED_LIBS=YES \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DENABLE_FORTRAN=ON \
        -DENABLE_PYTHON=OFF \
        -DBUILD_TESTING=OFF .. \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ../..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi
    fi
    LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXC from system paths ===================="
    check_lib -lxcf03 "libxc"
    check_lib -lxc "libxc"
    add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
    add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking LIBXC to user paths ===================="
    pkg_install_dir="$with_libxc"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libxc" != "__DONTUSE__" ]; then
  LIBXC_LIBS="-lxcf03 -lxc"
  if [ "$with_libxc" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libxc"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export LD_LIBRARY_PATH="$pkg_install_dir/lib":\${LD_LIBRARY_PATH}
export LD_RUN_PATH="$pkg_install_dir/lib":\${LD_RUN_PATH}
export LIBRARY_PATH="$pkg_install_dir/lib":\${LIBRARY_PATH}
export CPATH="$pkg_install_dir/include":\${CPATH}
export PKG_CONFIG_PATH="$pkg_install_dir/lib/pkgconfig":\${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH="$pkg_install_dir":\${CMAKE_PREFIX_PATH}
EOF
    cat "${BUILDDIR}/setup_libxc" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libxc"
export LIBXC_CFLAGS="${LIBXC_CFLAGS}"
export LIBXC_LDFLAGS="${LIBXC_LDFLAGS}"
export LIBXC_LIBS="${LIBXC_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBXC"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBXC_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBXC_LDFLAGS}"
export CP_LIBS="${LIBXC_LIBS} \${CP_LIBS}"
export LIBXC_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_libxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libxc"
