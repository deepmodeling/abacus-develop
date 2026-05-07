#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-01-04
# other contributor: Benrui Tang

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${SCRIPT_DIR}"/package_versions.sh

# Load ELPA package variables with version suffix support
# Check for version configuration from environment or individual package setting
version_suffix=""
if [[ -n "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" ]]; then
    # Check for individual package version override
    if echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "elpa:alt"; then
        version_suffix="alt"
    elif echo "${ABACUS_TOOLCHAIN_PACKAGE_VERSIONS}" | grep -q "elpa:main"; then
        version_suffix="main"
    fi
fi
# Fall back to global version suffix if no individual setting
if [[ -z "$version_suffix" && -n "${ABACUS_TOOLCHAIN_VERSION_SUFFIX}" ]]; then
    version_suffix="${ABACUS_TOOLCHAIN_VERSION_SUFFIX}"
fi
# Load package variables with appropriate version
load_package_vars "elpa" "$version_suffix"
elpa_pkg="elpa-${elpa_ver}.tar.gz"

source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_elpa" ] && rm "${BUILDDIR}/setup_elpa"

ELPA_CFLAGS=''
ELPA_LDFLAGS=''
ELPA_LIBS=''
elpa_dir_openmp="_openmp"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# elpa only works with MPI switched on
if [ $MPI_MODE = no ]; then
    report_warning $LINENO "MPI is disabled, skipping elpa installation"
    exit 0
fi

case "$with_elpa" in
    __INSTALL__)
        echo "==================== Installing ELPA ===================="
        pkg_install_dir="${INSTALLDIR}/elpa-${elpa_ver}"
        #pkg_install_dir="${HOME}/lib/elpa/${elpa_ver}-gcc8"
        install_lock_file="$pkg_install_dir/install_successful"
        enable_openmp="yes"

        # specific settings needed on CRAY Linux Environment
        if [ "$ENABLE_CRAY" = "__TRUE__" ]; then
            if [ ${CRAY_PRGENVCRAY} ]; then
                # extra LDFLAGS needed
                cray_ldflags="-dynamic"
            fi
            # enable_openmp="no"
        fi

        if verify_checksums "${install_lock_file}"; then
            echo "elpa-${elpa_ver} is already installed, skipping it."
        else
            require_env MATH_LIBS
            url="https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${elpa_ver}/${elpa_pkg}"
            if [ -f ${elpa_pkg} ]; then
                echo "${elpa_pkg} is found"
            else
                download_pkg_from_url "${elpa_sha256}" "${elpa_pkg}" "${url}"
            fi
            if [ "${PACK_RUN}" = "__TRUE__" ]; then
                echo "--pack-run mode specified, skip installation"
                exit 0
            fi
            [ -d elpa-${elpa_ver} ] && rm -rf elpa-${elpa_ver}
            tar -xzf ${elpa_pkg}

            # elpa expect FC to be an mpi fortran compiler that is happy
            # with long lines, and that a bunch of libs can be found
            cd elpa-${elpa_ver}

            # ELPA-2017xxxx enables AVX2 by default, switch off if machine doesn't support it.
            AVX_flag=""
            AVX512_flags=""
            FMA_flag=""
            SSE4_flag=""
            config_flags="--disable-avx-kernels --disable-avx2-kernels --disable-avx512-kernels --disable-sse-kernels --disable-sse-assembly-kernels"
            if [ "${TARGET_CPU}" = "native" ]; then
                if [ -f /proc/cpuinfo ] && [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
                    has_AVX=$(grep '\bavx\b' /proc/cpuinfo 1> /dev/null && echo 'yes' || echo 'no')
                    [ "${has_AVX}" = "yes" ] && AVX_flag="-mavx" || AVX_flag=""
                    has_AVX2=$(grep '\bavx2\b' /proc/cpuinfo 1> /dev/null && echo 'yes' || echo 'no')
                    [ "${has_AVX2}" = "yes" ] && AVX_flag="-mavx2"
                    has_AVX512=$(grep '\bavx512f\b' /proc/cpuinfo 1> /dev/null && echo 'yes' || echo 'no')
                    [ "${has_AVX512}" = "yes" ] && AVX512_flags="-mavx512f"
                    FMA_flag=$(grep '\bfma\b' /proc/cpuinfo 1> /dev/null && echo '-mfma' || echo '-mno-fma')
                    SSE4_flag=$(grep '\bsse4_1\b' /proc/cpuinfo 1> /dev/null && echo '-msse4' || echo '-mno-sse4')
                    grep '\bavx512dq\b' /proc/cpuinfo 1> /dev/null && AVX512_flags+=" -mavx512dq"
                    grep '\bavx512cd\b' /proc/cpuinfo 1> /dev/null && AVX512_flags+=" -mavx512cd"
                    grep '\bavx512bw\b' /proc/cpuinfo 1> /dev/null && AVX512_flags+=" -mavx512bw"
                    grep '\bavx512vl\b' /proc/cpuinfo 1> /dev/null && AVX512_flags+=" -mavx512vl"
                    config_flags="--enable-avx-kernels=${has_AVX} --enable-avx2-kernels=${has_AVX2} --enable-avx512-kernels=${has_AVX512}"
                fi
            fi
            for TARGET in "cpu" "nvidia"; do
                # Accept both uppercase and lowercase GPU enable flags for compatibility
                gpu_enabled="${ENABLE_CUDA:-${enable_cuda}}"
                [ "$TARGET" = "nvidia" ] && [ "$gpu_enabled" != "__TRUE__" ] && continue
                # disable cpu if cuda is enabled, only install one
                [ "$TARGET" != "nvidia" ] && [ "$gpu_enabled" = "__TRUE__" ] && continue
                # extend the pkg_install_dir by TARGET
                # this linking method is totally different from cp2k toolchain
                # for cp2k, ref https://github.com/cp2k/cp2k/commit/6fe2fc105b8cded84256248f68c74139dd8fc2e9
                pkg_install_dir="${pkg_install_dir}/${TARGET}"

                echo "Installing from scratch into ${pkg_install_dir}"
                mkdir -p "build_${TARGET}"
                cd "build_${TARGET}"
                if [ "${with_amd}" != "__DONTUSE__" ] && [ "${WITH_FLANG}" = "yes" ] ; then
                    # special option for flang compiler
                    echo "AMD fortran compiler detected, enable special option operation"
                    ../configure --prefix="${pkg_install_dir}" \
                        --libdir="${pkg_install_dir}/lib" \
                        --enable-openmp=${enable_openmp} \
                        --enable-static=no \
                        --enable-shared=yes \
                        --disable-c-tests \
                        --disable-cpp-tests \
                        ${config_flags} \
                        --with-cuda-path=${CUDA_PATH:-${CUDA_HOME:-/CUDA_HOME-notset}} \
                        --enable-nvidia-gpu-kernels=$([ "$TARGET" = "nvidia" ] && echo "yes" || echo "no") \
                        --with-NVIDIA-GPU-compute-capability=$([ "$TARGET" = "nvidia" ] && echo "sm_$ARCH_NUM" || echo "sm_70") \
                        --enable-nvidia-cub --with-cusolver \
                        OMPI_MCA_plm_rsh_agent=/bin/false \
                        FC=${MPIFC} \
                        CC=${MPICC} \
                        CXX=${MPICXX} \
                        CPP="cpp -E" \
                        FCFLAGS="${FCFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        CFLAGS="${CFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        CXXFLAGS="${CXXFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        LDFLAGS="${MATH_LDFLAGS} ${SCALAPACK_LDFLAGS} ${cray_ldflags} -lstdc++" \
                        LIBS="${SCALAPACK_LIBS} $(resolve_string "${MATH_LIBS}" "MPI") ${MPI_LIBS}" \
                        SCALAPACK_LDFLAGS="${SCALAPACK_LDFLAGS}" \
                        SCALAPACK_FCFLAGS="${SCALAPACK_CFLAGS}" \
                        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
                    # remove unsupported compile option in libtool
                    sed -i ./libtool \
                        -e 's/\\$wl-soname //g' \
                        -e 's/\\$wl--whole-archive\\$convenience \\$wl--no-whole-archive//g' \
                        -e 's/\\$wl\\$soname //g'
                else
                    # normal installation
                    ../configure --prefix="${pkg_install_dir}/" \
                        --libdir="${pkg_install_dir}/lib" \
                        --enable-openmp=${enable_openmp} \
                        --enable-static=no \
                        --enable-shared=yes \
                        --disable-c-tests \
                        --disable-cpp-tests \
                        ${config_flags} \
                        --enable-nvidia-gpu-kernels=$([ "$TARGET" = "nvidia" ] && echo "yes" || echo "no") \
                        --with-cuda-path=${CUDA_PATH:-${CUDA_HOME:-/CUDA_HOME-notset}} \
                        --with-NVIDIA-GPU-compute-capability=$([ "$TARGET" = "nvidia" ] && echo "sm_$ARCH_NUM" || echo "sm_70") \
                        --enable-nvidia-cub --with-cusolver \
                        FC=${MPIFC} \
                        CC=${MPICC} \
                        CXX=${MPICXX} \
                        CPP="cpp -E" \
                        FCFLAGS="${FCFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        CFLAGS="${CFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        CXXFLAGS="${CXXFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} ${AVX_flag} ${FMA_flag} ${SSE4_flag} ${AVX512_flags} -fno-lto" \
                        LDFLAGS="-Wl,--allow-multiple-definition -Wl,--enable-new-dtags ${MATH_LDFLAGS} ${SCALAPACK_LDFLAGS} ${cray_ldflags} -lstdc++" \
                        LIBS="${SCALAPACK_LIBS} $(resolve_string "${MATH_LIBS}" "MPI")" \
                        SCALAPACK_LDFLAGS="${SCALAPACK_LDFLAGS}" \
                        SCALAPACK_FCFLAGS="${SCALAPACK_CFLAGS}" \
                        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
                fi
                make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
                make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
                cd ..
                # link elpa
                link=${pkg_install_dir}/include/elpa
                if [[ ! -d $link ]]; then
                    ln -s ${pkg_install_dir}/include/elpa_openmp-${elpa_ver}/elpa $link
                fi
            done
            cd ..
            
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
        fi
        [ "$enable_openmp" != "yes" ] && elpa_dir_openmp=""
        ELPA_CFLAGS="-I'${pkg_install_dir}/include/elpa${elpa_dir_openmp}-${elpa_ver}/modules' -I'${pkg_install_dir}/include/elpa${elpa_dir_openmp}-${elpa_ver}/elpa'"
        ELPA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding ELPA from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        check_lib -lelpa_openmp "ELPA"
        # get the include paths
        elpa_include="$(find_in_paths "elpa_openmp-*" $INCLUDE_PATHS)"
        if [ "$elpa_include" != "__FALSE__" ]; then
            echo "ELPA include directory threaded version is found to be $elpa_include"
            ELPA_CFLAGS="-I'$elpa_include/modules' -I'$elpa_include/elpa'"
        else
            echo "Cannot find elpa_openmp-${elpa_ver} from paths $INCLUDE_PATHS"
            exit 1
        fi
        # get the lib paths
        add_lib_from_paths ELPA_LDFLAGS "libelpa.*" $LIB_PATHS
        ;;
    __DONTUSE__) ;;

    *)
        echo "==================== Linking ELPA to user paths ===================="
        pkg_install_dir="$with_elpa"
        check_dir "${pkg_install_dir}/include"
        check_dir "${pkg_install_dir}/lib"
        user_include_path="${pkg_install_dir}/include"
        elpa_include="$(find_in_paths "elpa_openmp-*" user_include_path)"
        if [ "$elpa_include" != "__FALSE__" ]; then
            echo "ELPA include directory threaded version is found to be $elpa_include/modules"
            check_dir "$elpa_include/modules"
            ELPA_CFLAGS="-I'$elpa_include/modules' -I'$elpa_include/elpa'"
        else
            echo "Cannot find elpa_openmp-* from path $user_include_path"
            exit 1
        fi
        ELPA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_elpa" != "__DONTUSE__" ]; then
    ELPA_LIBS="-lelpa${elpa_dir_openmp}"
    cat << EOF > "${BUILDDIR}/setup_elpa"
prepend_path CPATH "$elpa_include"
EOF
    if [ "$with_elpa" != "__SYSTEM__" ]; then
        cat << EOF >> "${BUILDDIR}/setup_elpa"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export PATH="$pkg_install_dir/bin":\${PATH}
export LD_LIBRARY_PATH="$pkg_install_dir/lib":\${LD_LIBRARY_PATH}
export LD_RUN_PATH="$pkg_install_dir/lib":\${LD_RUN_PATH}
export LIBRARY_PATH="$pkg_install_dir/lib":\${LIBRARY_PATH}
export CPATH="$pkg_install_dir/include":\${CPATH}
export PKG_CONFIG_PATH="$pkg_install_dir/lib/pkgconfig":\${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH="$pkg_install_dir":\${CMAKE_PREFIX_PATH}
export ELPA_ROOT="$pkg_install_dir"
EOF
        cat "${BUILDDIR}/setup_elpa" >> $SETUPFILE
    fi
    cat << EOF >> "${BUILDDIR}/setup_elpa"
export ELPA_CFLAGS="${ELPA_CFLAGS}"
export ELPA_LDFLAGS="${ELPA_LDFLAGS}"
export ELPA_LIBS="${ELPA_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__ELPA IF_CUDA(-D__ELPA_NVIDIA_GPU|)|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${ELPA_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${ELPA_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${ELPA_LIBS}|) \${CP_LIBS}"
export ELPA_ROOT="$pkg_install_dir"
export ELPA_VERSION="${elpa_ver}"
EOF

fi

load "${BUILDDIR}/setup_elpa"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "elpa"
