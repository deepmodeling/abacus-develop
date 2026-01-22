#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-0504

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load CMake package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "cmake:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "cmake:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi

# Ensure OPENBLAS_ARCH is set before loading package variables
# This is needed for architecture-specific SHA256 selection
# In --pack-run mode, openblas_arch file may contain empty values, so we need fallback
if [ -f "${BUILDDIR}/openblas_arch" ]; then
    source "${BUILDDIR}/openblas_arch"
fi

if [ -z "${OPENBLAS_ARCH}" ]; then
    case "$(uname -m)" in
        x86_64|amd64) OPENBLAS_ARCH="x86_64" ;;
        aarch64|arm64) OPENBLAS_ARCH="arm64" ;;
        *) OPENBLAS_ARCH="x86_64" ;;  # default fallback
    esac
    echo "OPENBLAS_ARCH not set, using fallback: ${OPENBLAS_ARCH}"
fi

# Export OPENBLAS_ARCH to ensure it's available throughout the script
export OPENBLAS_ARCH

# Load package variables with appropriate version
load_package_vars "cmake" "$version_suffix"

source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# Re-apply architecture detection if OPENBLAS_ARCH is still empty after sourcing
if [ -z "${OPENBLAS_ARCH}" ]; then
    case "$(uname -m)" in
        x86_64|amd64) OPENBLAS_ARCH="x86_64" ;;
        aarch64|arm64) OPENBLAS_ARCH="arm64" ;;
        *) OPENBLAS_ARCH="x86_64" ;;  # default fallback
    esac
    export OPENBLAS_ARCH
fi

[ -f "${BUILDDIR}/setup_cmake" ] && rm "${BUILDDIR}/setup_cmake"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "${with_cmake}" in
    __INSTALL__)
        echo "==================== Installing CMake ===================="
        cmake_ext="sh"
        if [ "${OPENBLAS_ARCH}" = "arm64" ]; then
            if [ "$(uname -s)" = "Darwin" ]; then
                cmake_arch="macos-universal"
                cmake_ext="tar.gz"
            elif [ "$(uname -s)" = "Linux" ]; then
                cmake_arch="linux-aarch64"
            else
                report_error ${LINENO} \
                    "cmake installation for ARCH=${OPENBLAS_ARCH} under $(uname -s) is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
            fi
        elif [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
            cmake_arch="linux-x86_64"
        else
            report_error ${LINENO} \
                "cmake installation for ARCH=${OPENBLAS_ARCH} is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
            exit 1
        fi
        pkg_install_dir="${INSTALLDIR}/cmake-${cmake_ver}"
        #pkg_install_dir="${HOME}/apps/cmake/${cmake_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        cmake_pkg="cmake-${cmake_ver}-${cmake_arch}.${cmake_ext}"
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
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            mkdir -p ${pkg_install_dir}
            if [ "${cmake_arch}" = "macos-universal" ]; then
                tar --strip-components=3 -xvf $cmake_pkg -C ${pkg_install_dir} > install.log 2>&1 || tail -n ${LOG_LINES} install.log
            else
                /bin/sh $cmake_pkg --prefix=${pkg_install_dir} --skip-license > install.log 2>&1 || tail -n ${LOG_LINES} install.log
            fi
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage0/$(basename ${SCRIPT_NAME})"
        fi
        ;;
    __SYSTEM__)
        echo "==================== Finding CMake from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
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
