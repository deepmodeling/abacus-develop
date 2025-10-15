#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_mkl" ] && rm "${BUILDDIR}/setup_mkl"

MKL_CFLAGS=""
MKL_LDFLAGS=""
MKL_LIBS=""
MKL_FFTW="yes"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_mkl}" in
    __INSTALL__)
        echo "==================== Installing MKL ===================="
        report_error ${LINENO} "To install MKL, please contact your system administrator."
        exit 1
        ;;
    __SYSTEM__)
        echo "==================== Finding MKL from system paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        if ! [ -z "${MKLROOT}" ]; then
            echo "MKLROOT is found to be ${MKLROOT}"
        else
            report_error ${LINENO} "Cannot find env variable MKLROOT, the script relies on it being set. Please check in MKL installation and use --with-mkl=<location> to pass the path to MKL root directory to this script."
            exit 1
        fi
        check_lib -lm
        check_lib -ldl
        ;;
    __DONTUSE__)
        # Nothing to do
        ;;
    *)
        echo "==================== Linking MKL to user paths ===================="
        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip system check"
            exit 0
        fi
        check_dir "${with_mkl}"
        MKLROOT="${with_mkl}"
        ;;
esac
if [ "${with_mkl}" != "__DONTUSE__" ]; then
    case ${OPENBLAS_ARCH} in
        x86_64)
            mkl_arch_dir="intel64"
            MKL_CFLAGS="-m64"
            ;;
        i386)
            mkl_arch_dir="ia32"
            MKL_CFLAGS="-m32"
            ;;
        *)
            report_error $LINENO "MKL only supports intel64 (x86_64) and ia32 (i386) at the moment, and your arch obtained from OpenBLAS prebuild is $OPENBLAS_ARCH"
            exit 1
            ;;
    esac
    mkl_lib_dir="${MKLROOT}/lib/${mkl_arch_dir}"
    # check we have required libraries
    mkl_required_libs="libmkl_gf_lp64.so libmkl_sequential.so libmkl_core.so"
    for ii in $mkl_required_libs; do
        if [ ! -f "$mkl_lib_dir/${ii}" ]; then
            report_error $LINENO "missing MKL library ${ii}"
            exit 1
        fi
    done
    MKL_CFLAGS="${MKL_CFLAGS} -I'${MKLROOT}/include'"
    MKL_LDFLAGS="-L'${mkl_lib_dir}' -Wl,-rpath,'${mkl_lib_dir}'"
    MKL_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
    cat << EOF > "${BUILDDIR}/setup_mkl"
prepend_path LD_LIBRARY_PATH "${mkl_lib_dir}"
prepend_path LD_RUN_PATH "${mkl_lib_dir}"
prepend_path LIBRARY_PATH "${mkl_lib_dir}"
prepend_path CPATH "${MKLROOT}/include"
export MKL_CFLAGS="${MKL_CFLAGS}"
export MKL_LDFLAGS="${MKL_LDFLAGS}"
export MKL_LIBS="${MKL_LIBS}"
export FAST_MATH_CFLAGS="\${FAST_MATH_CFLAGS} ${MKL_CFLAGS}"
export FAST_MATH_LDFLAGS="\${FAST_MATH_LDFLAGS} ${MKL_LDFLAGS}"
export FAST_MATH_LIBS="\${FAST_MATH_LIBS} ${MKL_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_BLAS(-D__MKL|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_BLAS(${MKL_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_BLAS(${MKL_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_BLAS(${MKL_LIBS}|)"
export MKLROOT="${MKLROOT}"
EOF
    cat "${BUILDDIR}/setup_mkl" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_mkl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "mkl"
