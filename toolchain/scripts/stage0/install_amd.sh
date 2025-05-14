#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-0308
# NOTICE: flang cannot be used when compiling ELPA

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=${0}
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_amd" ] && rm "${BUILDDIR}/setup_amd"

AMD_CFLAGS=""
AMD_LDFLAGS=""
AMD_LIBS=""
mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_amd}" in
  __INSTALL__)
    echo "==================== Installing the AMD compiler ======================"
    echo "__INSTALL__ is not supported; please install the AMD compiler manually"
    if [ "${PACK_RUN}" != "__TRUE__" ]; then
        exit 1
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding AMD compiler from system paths ===================="
    if [ "${PACK_RUN}" = "__TRUE__" ]; then
        echo "--pack-run mode specified, skip system check"
    else
    check_command clang "amd" && CC="$(realpath $(command -v clang))" || exit 1
    check_command clang++ "amd" && CXX="$(realpath $(command -v clang++))" || exit 1
    if [ "${WITH_FLANG}" = "yes" ]; then
        check_command flang "amd" && FC="$(realpath $(command -v flang))" || exit 1
    else
        check_command gfortran "gcc" && FC="gfortran" || exit 1
        add_lib_from_paths GCC_LDFLAGS "libgfortran.*" ${LIB_PATHS}
    fi
    fi
    F90="${FC}"
    F77="${FC}"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking AMD compiler to user paths ===================="
    pkg_install_dir="${with_amd}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/clang "amd" && CC="${pkg_install_dir}/bin/clang" || exit 1
    check_command ${pkg_install_dir}/bin/clang++ "amd" && CXX="${pkg_install_dir}/bin/clang++" || exit 1
    if [ "${WITH_FLANG}" = "yes" ]; then
        check_command ${pkg_install_dir}/bin/flang "amd" && FC="${pkg_install_dir}/bin/flang" || exit 1
    else
        check_command gfortran "gcc" && FC="$(command -v gfortran)" || exit 1
    fi
    F90="${FC}"
    F77="${FC}"
    AMD_CFLAGS="-I'${pkg_install_dir}/include'"
    AMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "${with_amd}" != "__DONTUSE__" ]; then
  echo "CC  is ${CC}"
  [ $(realpath $(command -v clang) | grep -e aocc-compiler) ] || echo "Check the AMD C compiler path"
  echo "CXX is ${CXX}"
  [ $(realpath $(command -v clang++) | grep -e aocc-compiler) ] || echo "Check the AMD C++ compiler path"
  echo "FC  is ${FC}"
  if [ "${WITH_FLANG}" = "yes" ]; then
    [ $(realpath $(command -v flang) | grep -e aocc-compiler) ] || echo "Check the AMD Fortran compiler path"
  else
    [ $(realpath $(command -v gfortran) | grep -e aocc-compiler) ] || echo "Check the GNU Fortran compiler path"
  fi
  cat << EOF > "${BUILDDIR}/setup_amd"
export CC="${CC}"
export CXX="${CXX}"
export FC="${FC}"
export F90="${F90}"
export F77="${F77}"
EOF
  if [ "${with_amd}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_amd"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_amd"
export AMD_CFLAGS="${AMD_CFLAGS}"
export AMD_LDFLAGS="${AMD_LDFLAGS}"
export AMD_LIBS="${AMD_LIBS}"
EOF
  cat "${BUILDDIR}/setup_amd" >> ${SETUPFILE}
fi

load "${BUILDDIR}/setup_amd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "amd"
