#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-10-10

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

nep_ver="main"
nep_sha256="--no-checksum"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_nep" ] && rm "${BUILDDIR}/setup_nep"

NEP_CFLAGS=""
NEP_LDFLAGS=""
NEP_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_nep" in
  __INSTALL__)
    echo "==================== Installing NEP (CPU version) ===================="
    dirname="NEP_CPU-${nep_ver}"
    pkg_install_dir="${INSTALLDIR}/${dirname}"
    install_lock_file="$pkg_install_dir/install_successful"
    filename="nep-${nep_ver}.tar.gz"
    url="https://codeload.github.com/brucefan1983/NEP_CPU/tar.gz/${nep_ver}"

    if verify_checksums "${install_lock_file}"; then
        echo "$dirname is already installed, skipping it."
    else
        if [ -f $filename ]; then
            echo "$filename is found"
        else
            echo "===> Notice: This version of NEP_CPU is downloaded from the GitHub master repository <==="
            download_pkg_from_url "${nep_sha256}" "${filename}" "${url}"
        fi

        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip installation"
        else
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d $dirname ] && rm -rf $dirname
            tar -xzf $filename
            cd $dirname

            cat << EOF > Makefile
CXX ?= g++

# Compiler flags
CXXFLAGS = -O2 -fPIC -std=c++11

# Include directories
INCLUDES = -I./src

# Source files
SRCS = ./src/nep.cpp

# Object files
OBJS = \$(SRCS:.cpp=.o)

# Target shared library
TARGET = libnep.so

# Default target
all: \$(TARGET)

# Rule to build the shared library
\$(TARGET): \$(OBJS)
	\$(CXX) -shared \$(OBJS) -o \$(TARGET)

# Rule to compile source files into object files
%.o: %.cpp
	\$(CXX) \$(CXXFLAGS) \$(INCLUDES) -c \$< -o \$@

# Clean up build files
clean:
	rm -f \$(OBJS) \$(TARGET)

# Install target
install:
	mkdir -p \$(PREFIX)/lib
	mkdir -p \$(PREFIX)/include
	cp \$(TARGET) \$(PREFIX)/lib/
	cp src/nep.h \$(PREFIX)/include/
EOF

            make > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            make PREFIX="${pkg_install_dir}" install > install.log 2>&1 || tail -n ${LOG_LINES} install.log

            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
        fi
    fi
    NEP_CFLAGS="-I'${pkg_install_dir}/include'"
    NEP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;

  __SYSTEM__)
    echo "==================== Finding NEP_CPU from system paths ===================="
    check_lib -lnep "nep"
    add_include_from_paths NEP_CFLAGS "nep.h" $INCLUDE_PATHS
    add_lib_from_paths NEP_LDFLAGS "libnep.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking NEP_CPU to user paths ===================="
    pkg_install_dir="$with_nep"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    NEP_CFLAGS="-I'${pkg_install_dir}/include'"
    NEP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_nep" != "__DONTUSE__" ]; then
  NEP_LIBS="-lnep"
  if [ "$with_nep" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_nep"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export LD_LIBRARY_PATH="$pkg_install_dir/lib":\${LD_LIBRARY_PATH}
export LD_RUN_PATH="$pkg_install_dir/lib":\${LD_RUN_PATH}
export LIBRARY_PATH="$pkg_install_dir/lib":\${LIBRARY_PATH}
export CPATH="$pkg_install_dir/include":\${CPATH}
export CMAKE_PREFIX_PATH="$pkg_install_dir":\${CMAKE_PREFIX_PATH}
EOF
    cat "${BUILDDIR}/setup_nep" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_nep"
export NEP_CFLAGS="${NEP_CFLAGS}"
export NEP_LDFLAGS="${NEP_LDFLAGS}"
export NEP_LIBS="${NEP_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__NEP"
export CP_CFLAGS="\${CP_CFLAGS} \${NEP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} \${NEP_LDFLAGS}"
export CP_LIBS="\${NEP_LIBS} \${CP_LIBS}"
export NEP_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_nep"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "nep"