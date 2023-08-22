#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# CEREAL is not need any complex setting

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

cereal_ver="6.2.2"
#cereal_sha256="a0f6f1bba7ba5c0c85b2bfe65aca1591025f509a7f11471b4cd651a79491b045"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cereal" ] && rm "${BUILDDIR}/setup_cereal"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_cereal" in
  __INSTALL__)
    echo "==================== Installing CEREAL ===================="
    pkg_install_dir="${INSTALLDIR}/cereal-${cereal_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "cereal-${cereal_ver} is already installed, skipping it."
    else
      if [ -f cereal-${cereal_ver}.tar.gz ]; then
        echo "cereal-${cereal_ver}.tar.gz is found"
      else
        download_pkg_from_ABACUS_org "${cereal_sha256}" "cereal-${cereal_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d cereal-${cereal_ver} ] && rm -rf cereal-${cereal_ver}
      tar -xzf cereal-${cereal_ver}.tar.gz
      cd cereal-${cereal_ver}

      # ABACUS does not make use of fourth derivatives, so skip their compilation with --disable-lxc
      ./configure --prefix="${pkg_install_dir}" --libdir="${pkg_install_dir}/lib" --enable-shared --enable-static \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi
    cereal_CFLAGS="-I'${pkg_install_dir}/include'"
    cereal_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding cereal from system paths ===================="
    check_lib -lxcf03 "cereal"
    check_lib -lxc "cereal"
    add_include_from_paths cereal_CFLAGS "xc.h" $INCLUDE_PATHS
    add_lib_from_paths cereal_LDFLAGS "cereal.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking cereal to user paths ===================="
    pkg_install_dir="$with_cereal"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    cereal_CFLAGS="-I'${pkg_install_dir}/include'"
    cereal_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_cereal" != "__DONTUSE__" ]; then
  cereal_LIBS="-lxcf03 -lxc"
  if [ "$with_cereal" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_cereal"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    cat "${BUILDDIR}/setup_cereal" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_cereal"
export cereal_CFLAGS="${cereal_CFLAGS}"
export cereal_LDFLAGS="${cereal_LDFLAGS}"
export cereal_LIBS="${cereal_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__cereal"
export CP_CFLAGS="\${CP_CFLAGS} ${cereal_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${cereal_LDFLAGS}"
export CP_LIBS="${cereal_LIBS} \${CP_LIBS}"
export cereal_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_cereal"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cereal"
