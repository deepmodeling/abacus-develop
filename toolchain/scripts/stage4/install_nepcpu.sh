#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Last Update in 2025-10-10

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

nepcpu_ver="main"
nepcpu_sha256="--no-checksum"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_nepcpu" ] && rm "${BUILDDIR}/setup_nepcpu"

NEPCPU_CFLAGS=""
NEPCPU_LDFLAGS=""
NEPCPU_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_nepcpu" in
  __INSTALL__)
    echo "==================== Installing NEP_CPU ===================="
    dirname="NEP_CPU-${nepcpu_ver}"
    pkg_install_dir="${INSTALLDIR}/nep_cpu-${nepcpu_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    filename="NEP_CPU-${nepcpu_ver}.tar.gz"
    url="https://codeload.github.com/brucefan1983/NEP_CPU/tar.gz/${nepcpu_ver}"

    if verify_checksums "${install_lock_file}"; then
        echo "$dirname is already installed, skipping it."
    else
        if [ -f $filename ]; then
            echo "$filename is found"
        else
            echo "===> Notice: This version of NEP_CPU is downloaded from the GitHub master repository <==="
            download_pkg_from_url "${nepcpu_sha256}" "${filename}" "${url}"
        fi

        if [ "${PACK_RUN}" = "__TRUE__" ]; then
            echo "--pack-run mode specified, skip installation"
        else
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d $dirname ] && rm -rf $dirname
            tar -xzf $filename
            cd $dirname

            cat << EOF > Makefile
# Compiler - Use CXX from the environment
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
TARGET = libnepcpu.so

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
    NEPCPU_CFLAGS="-I'${pkg_install_dir}/include'"
    NEPCPU_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;

  __SYSTEM__)
    echo "==================== Finding NEP_CPU from system paths ===================="
    echo "Finding NEP_CPU from system paths is not supported yet."
    echo "Please use __INSTALL__ or provide a path to a pre-installed version."
    exit 1
    ;;

  __DONTUSE__) ;;

  *)
    echo "==================== Linking NEP_CPU to user paths ===================="
    pkg_install_dir="$with_nepcpu"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    NEPCPU_CFLAGS="-I'${pkg_install_dir}/include'"
    NEPCPU_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_nepcpu" != "__DONTUSE__" ]; then
  NEPCPU_LIBS="-lnepcpu"
  if [ "$with_nepcpu" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_nepcpu"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
export LD_LIBRARY_PATH="$pkg_install_dir/lib":\${LD_LIBRARY_PATH}
export LD_RUN_PATH="$pkg_install_dir/lib":\${LD_RUN_PATH}
export LIBRARY_PATH="$pkg_install_dir/lib":\${LIBRARY_PATH}
export CPATH="$pkg_install_dir/include":\${CPATH}
EOF
    cat "${BUILDDIR}/setup_nepcpu" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_nepcpu"
export NEPCPU_CFLAGS="${NEPCPU_CFLAGS}"
export NEPCPU_LDFLAGS="${NEPCPU_LDFLAGS}"
export NEPCPU_LIBS="${NEPCPU_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__NEPCPU"
export CP_CFLAGS="\${CP_CFLAGS} \${NEPCPU_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} \${NEPCPU_LDFLAGS}"
export CP_LIBS="\${NEPCPU_LIBS} \${CP_LIBS}"
export NEPCPU_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_nepcpu"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "nepcpu"