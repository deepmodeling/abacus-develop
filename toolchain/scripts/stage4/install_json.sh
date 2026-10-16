#!/bin/bash -e

SCRIPT_NAME="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

source "${SCRIPT_DIR}/common_vars.sh"
source "${SCRIPT_DIR}/tool_kit.sh"
source "${SCRIPT_DIR}/signal_trap.sh"
source "${SCRIPT_DIR}/package_versions.sh"
load_package_vars "json"
source "${INSTALLDIR}/toolchain.conf"
source "${INSTALLDIR}/toolchain.env"

rm -f "${BUILDDIR}/setup_json"
mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
pkg_install_dir=""

case "${with_json}" in
    __INSTALL__)
        echo "==================== Installing nlohmann-json ===================="
        filename="json-${json_ver}.tar.xz"
        url="https://github.com/nlohmann/json/releases/download/v${json_ver}/json.tar.xz"
        pkg_install_dir="${INSTALLDIR}/json-${json_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"
        if verify_checksums "${install_lock_file}"; then
            echo "json-${json_ver} is already installed, skipping it."
        else
            retrieve_package "${json_sha256}" "${filename}" "${url}"
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            [ -d "json-${json_ver}" ] && rm -rf "json-${json_ver}"
            tar -xJf "json-${json_ver}.tar.xz"
            cd "json"
            mkdir build && cd build
            cmake .. \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
              -DJSON_BuildTests=OFF \
              > cmake.log 2>&1 || tail -n "${LOG_LINES}" cmake.log
            make install -j "$(get_nprocs)" > install.log 2>&1 || tail -n "${LOG_LINES}" install.log
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/install_json.sh"
        fi
        ;;
    __SYSTEM__)
        # ABACUS resolves the installed configuration with CMake's find_package.
        echo "Using system nlohmann-json (resolved by CMake)."
        ;;
    __DONTUSE__) ;;
    *)
        pkg_install_dir="${with_json}"
        check_dir "${pkg_install_dir}"
        ;;
esac

if [ -n "${pkg_install_dir}" ]; then
    cat > "${BUILDDIR}/setup_json" <<EOF
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
    filter_setup "${BUILDDIR}/setup_json" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_json"
write_toolchain_env "${INSTALLDIR}"
cd "${ROOTDIR}"
report_timing "json"
