#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-0308

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_aocl" ] && rm "${BUILDDIR}/setup_aocl"

AOCL_CFLAGS=""
AOCL_LDFLAGS=""
AOCL_LIBS=""
AOCL_ROOT=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_aocl}" in
  __INSTALL__)
    echo "==================== Installing AOCL ===================="
    report_error ${LINENO} "To install AOCL, please contact your system administrator."
    if [ "${PACK_RUN}" != "__TRUE__" ]; then
        exit 1
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding AOCL from system paths ===================="
    if [ "${PACK_RUN}" = "__TRUE__" ]; then
        echo "--pack-run mode specified, skip system check"
    else
        check_lib -lblis "AOCL"
        check_lib -lflame "AOCL"
        AOCL_LIBS="-lblis -lflame"
        add_include_from_paths AOCL_CFLAGS "blis.h" $INCLUDE_PATHS
        add_lib_from_paths AOCL_LDFLAGS "libblis.*" $LIB_PATHS
        add_include_from_paths AOCL_CFLAGS "lapack.h" $INCLUDE_PATHS
        add_lib_from_paths AOCL_LDFLAGS "libflame.*" $LIB_PATHS
    fi
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking AOCL to user paths ===================="
    pkg_install_dir="$with_openblas"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    AOCL_CFLAGS="-I'${pkg_install_dir}/include'"
    AOCL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    AOCL_LIBS="-lblis -lflame"
    ;;
esac
if [ "$with_openblas" != "__DONTUSE__" ]; then
  if [ "$with_openblas" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_aocl"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
prepend_path CPATH "$pkg_install_dir/include"
export LD_LIBRARY_PATH="$pkg_install_dir/lib:"\${LD_LIBRARY_PATH}
export LD_RUN_PATH="$pkg_install_dir/lib:"\${LD_RUN_PATH}
export LIBRARY_PATH="$pkg_install_dir/lib:"\${LIBRARY_PATH}
export CPATH="$pkg_install_dir/include:"\${CPATH}
export PKG_CONFIG_PATH="$pkg_install_dir/lib/pkgconfig:"\${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH="$pkg_install_dir:"\${CMAKE_PREFIX_PATH}
export AOCL_ROOT=${pkg_install_dir}
EOF
    cat "${BUILDDIR}/setup_aocl" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_aocl"
export AOCL_ROOT="${pkg_install_dir}"
export AOCL_CFLAGS="${AOCL_CFLAGS}"
export AOCL_LDFLAGS="${AOCL_LDFLAGS}"
export AOCL_LIBS="${AOCL_LIBS}"
export MATH_CFLAGS="\${MATH_CFLAGS} ${AOCL_CFLAGS}"
export MATH_LDFLAGS="\${MATH_LDFLAGS} ${AOCL_LDFLAGS}"
export MATH_LIBS="\${MATH_LIBS} ${AOCL_LIBS}"
export PKG_CONFIG_PATH="${pkg_install_dir}/lib/pkgconfig"
export CMAKE_PREFIX_PATH="${pkg_install_dir}"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_aocl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "aocl"
