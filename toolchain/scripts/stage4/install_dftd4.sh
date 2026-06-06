#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load package variables with appropriate version
load_package_vars "dftd4"
dftd4_pkg="dftd4-${dftd4_ver}.tar.xz"
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_dftd4" ] && rm "${BUILDDIR}/setup_dftd4"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_dftd4}" in
    __INSTALL__)
        echo "==================== Installing DFT-D4 ===================="
        url="https://github.com/dftd4/dftd4/releases/download/v${dftd4_ver}/${dftd4_pkg}"
        pkg_install_dir="${INSTALLDIR}/dftd4-${dftd4_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"

        if verify_checksums "${install_lock_file}"; then
            echo "dftd4-${dftd4_ver} is already installed, skipping it."
        else
            retrieve_package "${dftd4_sha256}" "${dftd4_pkg}" "${url}"
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi

            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d "dftd4-${dftd4_ver}" ] && rm -rf "dftd4-${dftd4_ver}"
            tar -xJf "${dftd4_pkg}"
            cd dftd4-${dftd4_ver}

            mkdir build && cd build
            CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${OPENBLAS_ROOT}" cmake \
                -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
                -DCMAKE_INSTALL_LIBDIR=lib \
                -DCMAKE_BUILD_TYPE=Release \
                -DCMAKE_VERBOSE_MAKEFILE=ON \
                .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
            make install -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
        ;;

    __SYSTEM__)
        echo "==================== Finding DFT-D4 from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        check_command pkg-config --modversion dftd4
        pkg_install_dir="$(dirname $(dirname $(find_in_paths "libdftd4.*" $LIB_PATHS)))"
        if [ -d "${pkg_install_dir}/lib/cmake/dftd4" ] ||
          [ -d "${pkg_install_dir}/lib64/cmake/dftd4" ]; then
          echo "Package dftd4 is found and confirmed to be built with CMake."
        else
          echo "ERROR: ABACUS requires dftd4 to be built with CMake."
          exit 1
        fi
        ;;

    __DONTUSE__) ;;

    *)
        echo "==================== Linking DFT-D4 to user paths ===================="
        pkg_install_dir="$with_dftd4"
        check_dir "${pkg_install_dir}/include"
        ;;
esac

if [ "${with_dftd4}" != "__DONTUSE__" ]; then
    cat << EOF > "${BUILDDIR}/setup_dftd4"
export DFTD4_VER="${dftd4_ver}"
export DFTD4_ROOT="${pkg_install_dir}"
EOF
    if [ "${with_dftd4}" != "__SYSTEM__" ]; then
        cat << EOF >> "${BUILDDIR}/setup_dftd4"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
    fi
    filter_setup "${BUILDDIR}/setup_dftd4" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_dftd4"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "dftd4"
