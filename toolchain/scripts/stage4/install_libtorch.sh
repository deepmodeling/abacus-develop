#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

# Last Update in 2024-0913

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# From https://pytorch.org/get-started/locally/
# libtorch_ver="1.12.1" 
# libtorch_sha256="82c7be80860f2aa7963f8700004a40af8205e1d721298f2e09b700e766a9d283"
# libtorch_ver="2.0.1" 
# libtorch_sha256="137a842d1cf1e9196b419390133a1623ef92f8f84dc7a072f95ada684f394afd"
libtorch_ver="2.1.2"
libtorch_sha256="904b764df6106a8a35bef64c4b55b8c1590ad9d071eb276e680cf42abafe79e9"

# user can manually download higher version of libtorch by:
# wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-{libtorch_ver}%2Bcpu.zip
# 2.4.0 latest, 2.1.2 recommended for lower GLIBC support (lower than 3.4.26)

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
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
        echo "Installing from scratch into ${pkg_install_dir}"
        [ -d libtorch ] && rm -rf libtorch
        [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
        unzip -q $filename # to libtorch
        mkdir -p "${pkg_install_dir}"
        mv libtorch/* "${pkg_install_dir}/"

        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    
    LIBTORCH_CXXFLAGS="-I${pkg_install_dir}/include"
    LIBTORCH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding libtorch from system paths ===================="
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
    else
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
  else
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
