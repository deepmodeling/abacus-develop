#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# contributor: MoseyQAQ (Denan Li)

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# Load version information from centralized package_versions.sh
source "${SCRIPT_DIR}/package_versions.sh"
load_package_vars "nep"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_nep" ] && rm "${BUILDDIR}/setup_nep"

NEP_CFLAGS=""
NEP_LDFLAGS=""
NEP_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_nep" in
  __INSTALL__)
    echo "==================== Installing NEP (CPU version) ===================="
    dirname="NEP_CPU-${nep_ver}"
    pkg_install_dir="${INSTALLDIR}/${dirname}"
    install_lock_file="${pkg_install_dir}/install_successful"
    filename="nep-${nep_ver}.tar.gz"
    url="https://codeload.github.com/brucefan1983/NEP_CPU/tar.gz/${nep_ver}"

    if verify_checksums "${install_lock_file}"; then
      echo "$dirname is already installed, skipping it."
    else
      retrieve_package "${nep_sha256}" "${filename}" "${url}"

      if [ "${PACK_RUN}" = "__TRUE__" ]; then
        echo "--pack-run mode specified, skip installation"
        exit 0
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d $dirname ] && rm -rf $dirname
      tar -xzf $filename
      cd $dirname
      cmake -Bbuild \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_BUILD_TYPE="Release" \
        -DNEP_CPU_OPENMP=ON \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      cmake --build build -t install -j "$(get_nprocs)" > build.log 2>&1 || tail -n ${LOG_LINES} build.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    NEP_CFLAGS="-I'${pkg_install_dir}/include'"
    NEP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;

  __SYSTEM__)
    echo "==================== Finding NEP_CPU from system paths ===================="
    check_lib -lnep "nep"
    add_include_from_paths NEP_CFLAGS "nep.h" $INCLUDE_PATHS
    add_lib_from_paths NEP_LDFLAGS "libnep.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking NEP_CPU to user paths ===================="
    pkg_install_dir="$with_nep"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    NEP_CFLAGS="-I'${pkg_install_dir}/include'"
    NEP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_nep" != "__DONTUSE__" ]; then
  NEP_LIBS="-lnep"
  if [ "$with_nep" != "__SYSTEM__" ]; then
    cat <<EOF >"${BUILDDIR}/setup_nep"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat <<EOF >>"${BUILDDIR}/setup_nep"
export NEP_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_nep" $SETUPFILE
fi

load "${BUILDDIR}/setup_nep"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "nep"
