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

# Load package variables with appropriate version
load_package_vars "cmake" "$version_suffix"

source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cmake" ] && rm "${BUILDDIR}/setup_cmake"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "${with_cmake}" in
    __INSTALL__)
        echo "==================== Installing CMake ===================="
        case "$(uname -s):${SYSTEM_ARCH}" in
            Darwin:x86_64 | Darwin:arm64)
                cmake_arch="macos-universal"
                cmake_checksum_arch="macos"
                ;;
            Linux:x86_64)
                cmake_arch="linux-x86_64"
                cmake_checksum_arch="x86_64"
                ;;
            Linux:arm64)
                cmake_arch="linux-aarch64"
                cmake_checksum_arch="aarch64"
                ;;
            *)
                report_error ${LINENO} \
                    "cmake installation for ARCH=${SYSTEM_ARCH} under $(uname -s) is not supported. You can try to use the system installation using the flag \"--with-cmake=system\" instead."
                exit 1
                ;;
        esac

        if [ "${version_suffix}" = "alt" ]; then
            cmake_checksum_var="cmake_alt_sha256_${cmake_checksum_arch}"
        else
            cmake_checksum_var="cmake_main_sha256_${cmake_checksum_arch}"
        fi
        cmake_sha256="${!cmake_checksum_var}"

        pkg_install_dir="${INSTALLDIR}/cmake-${cmake_ver}"
        #pkg_install_dir="${HOME}/apps/cmake/${cmake_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"
        cmake_pkg="cmake-${cmake_ver}-${cmake_arch}.tar.gz"
        if verify_checksums "${install_lock_file}"; then
            echo "cmake-${cmake_ver} is already installed, skipping it."
        else
            url="https://cmake.org/files/v${cmake_ver%.*}/${cmake_pkg}"
            retrieve_package "${cmake_sha256}" "${cmake_pkg}" "${url}"
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            mkdir -p ${pkg_install_dir}
            if [ "${cmake_arch}" = "macos-universal" ]; then
                strip_components=3
            else
                strip_components=1
            fi
            tar --strip-components=${strip_components} -xvf ${cmake_pkg} -C ${pkg_install_dir} > install.log 2>&1 || tail_excerpt install.log
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
EOF
        filter_setup "${BUILDDIR}/setup_cmake" $SETUPFILE
    fi
fi

load "${BUILDDIR}/setup_cmake"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cmake"
