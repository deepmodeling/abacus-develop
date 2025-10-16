#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

# Last Update in 2024-0913

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load LibTorch package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libtorch:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "libtorch:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "libtorch" "$version_suffix" 
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libtorch" ] && rm "${BUILDDIR}/setup_libtorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libtorch}" in
    __INSTALL__)
        echo "==================== Installing libtorch ===================="
        dirname="libtorch-${libtorch_ver}"
        pkg_install_dir="${INSTALLDIR}/${dirname}"
        #pkg_install_dir="${HOME}/lib/libtorch/${libtorch_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"
        archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}%2Bcpu.zip"
        filename="${dirname}.zip"

        if verify_checksums "${install_lock_file}"; then
            echo "${dirname} is already installed, skipping it."
        else
            if [ -f ${filename} ]; then
                echo "${filename} is found"
            else
                # download from pytorch.com and checksum
                url=https://download.pytorch.org/libtorch/cpu/${archive_file}
                download_pkg_from_url "${libtorch_sha256}" "${filename}" "${url}"
            fi
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            unzip -q $filename
            mkdir -p "${pkg_install_dir}"
            mv libtorch/* "${pkg_install_dir}/"

            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
        LIBTORCH_CXXFLAGS="-I${pkg_install_dir}/include"
        LIBTORCH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding libtorch from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        check_lib -ltorch "libtorch"
        add_include_from_paths LIBTORCH_CXXFLAGS "libtorch.h" $INCLUDE_PATHS
        add_lib_from_paths LIBTORCH_LDFLAGS "libtorch.*" "$LIB_PATHS"
        ;;
    __DONTUSE__) ;;

    *)
        echo "==================== Linking libtorch to user paths ===================="
        pkg_install_dir="${with_libtorch}"

        # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
        LIBTORCH_LIBDIR="${pkg_install_dir}/lib"
        [ -d "${pkg_install_dir}/lib64" ] && LIBTORCH_LIBDIR="${pkg_install_dir}/lib64"

        check_dir "${LIBTORCH_LIBDIR}"
        LIBTORCH_CXXFLAGS="-I${pkg_install_dir}/include"
        if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
            LIBTORCH_LDFLAGS="-Wl,--no-as-needed,-L'${LIBTORCH_LIBDIR}' -Wl,--no-as-needed,-rpath='${LIBTORCH_LIBDIR}'"
            LIBTORCH_LDFLAGS="-L'${LIBTORCH_LIBDIR}' -Wl,-rpath='${LIBTORCH_LIBDIR}'"
        fi
        ;;
esac

if [ "$with_libtorch" != "__DONTUSE__" ]; then
    if [ "$with_libtorch" != "__SYSTEM__" ]; then
        cat << EOF > "${BUILDDIR}/setup_libtorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export LD_LIBRARY_PATH="${pkg_install_dir}/lib":\${LD_LIBRARY_PATH}
export LD_RUN_PATH="${pkg_install_dir}/lib":\${LD_RUN_PATH}
export LIBRARY_PATH="${pkg_install_dir}/lib":\${LIBRARY_PATH}
export CPATH="${pkg_install_dir}/include":\${CPATH}
export PKG_CONFIG_PATH="${pkg_install_dir}/lib/pkgconfig":\${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH="${pkg_install_dir}":\${CMAKE_PREFIX_PATH}
EOF
    fi
    if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        cat << EOF >> "${BUILDDIR}/setup_libtorch"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBTORCH"
export CXXFLAGS="\${CXXFLAGS} ${LIBTORCH_CXXFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBTORCH_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} -lc10 -lc10_cuda -ltorch_cpu -ltorch_cuda -ltorch"
EOF
        cat "${BUILDDIR}/setup_libtorch" >> "${SETUPFILE}"
        cat << EOF >> "${BUILDDIR}/setup_libtorch"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBTORCH"
export CXXFLAGS="\${CXXFLAGS} ${LIBTORCH_CXXFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBTORCH_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} -lc10 -ltorch_cpu -ltorch"
EOF
        cat "${BUILDDIR}/setup_libtorch" >> "${SETUPFILE}"
    fi
fi

load "${BUILDDIR}/setup_libtorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libtorch"
