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

[ -f "${BUILDDIR}/setup_cmake" ] && rm "${BUILDDIR}/setup_cmake"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "${with_cmake}" in
  __INSTALL__)
    echo "==================== Installing CMake ===================="
    cmake_ver="3.31.2"
    if [ "${OPENBLAS_ARCH}" = "arm64" ]; then
      cmake_arch="linux-aarch64"
      cmake_sha256="85cc81f782cd8b5ac346e570ad5cfba3bdbe5aa01f27f7ce6266c4cef93342550"
    elif [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
      cmake_arch="linux-x86_64"
      cmake_sha256="b81cf3f4892683133f330cd7c016c28049b5725617db24ca8763360883545d34"
    else
      report_error ${LINENO} \
        "cmake installation for ARCH=${ARCH} is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
      exit 1
    fi
    pkg_install_dir="${INSTALLDIR}/cmake-${cmake_ver}"
    #pkg_install_dir="${HOME}/apps/cmake/${cmake_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    cmake_pkg="cmake-${cmake_ver}-${cmake_arch}.sh"
    if verify_checksums "${install_lock_file}"; then
      echo "cmake-${cmake_ver} is already installed, skipping it."
    else
      if [ -f $cmake_pkg ]; then
        echo "$cmake_pkg is found"
      else
        #download_pkg_from_ABACUS_org "${cmake_sha256}" "$cmake_pkg"
        url="https://cmake.org/files/v${cmake_ver%.*}/${cmake_pkg}"
        download_pkg_from_url "${cmake_sha256}" "${cmake_pkg}" "${url}"
      fi
      if [ "${PACK_RUN}" = "__TRUE__" ]; then
        echo "--pack-run mode specified, skip installation"
      else
        echo "Installing from scratch into ${pkg_install_dir}"
        mkdir -p ${pkg_install_dir}
        /bin/sh $cmake_pkg --prefix=${pkg_install_dir} --skip-license > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage0/$(basename ${SCRIPT_NAME})"
      fi
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding CMake from system paths ===================="
    check_command cmake "cmake"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking CMake to user paths ===================="
    pkg_install_dir="$with_cmake"
    check_dir "${with_cmake}/bin"
    ;;
esac
if [ "${with_cmake}" != "__DONTUSE__" ]; then
  if [ "${with_cmake}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_cmake"
prepend_path PATH "${pkg_install_dir}/bin"
export PATH="${pkg_install_dir}/bin":\${PATH}
EOF
    cat "${BUILDDIR}/setup_cmake" >> $SETUPFILE
  fi
fi

load "${BUILDDIR}/setup_cmake"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cmake"
