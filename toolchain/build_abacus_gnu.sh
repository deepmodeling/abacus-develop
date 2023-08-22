#!/bin/bash
# install ABACUS with libxc and deepks
# JamesMisaka in 2023.08.22

# Build ABACUS by gnu-toolchain

#rm -rf ../build
# module load openmpi
TOOL=$(pwd)
ABACUS_DIR=..
source ./install/setup
cd $ABACUS_DIR

PREFIX=.
BUILD_DIR=build_abacus
LAPACK=$TOOL/install/openblas-0.3.23/lib
SCALAPACK=$TOOL/install/scalapalack-2.2.1/lib
ELPA=$TOOL/install/elpa-2021.11.002/cpu
FFTW3=$TOOL/install/fftw-3.3.10
CEREAL=$TOOL/libs/cereal-1.3.2/
LIBXC=$TOOL/install/libxc-6.2.2
LIBTORCH=$TOOL/install/libtorch-2.0.1/share/cmake/Torch
LIBNPY=$TOOL/libs/libnpy-0.1.0/
DEEPMD=$HOME/apps/anaconda3/envs/deepmd

CC=gcc
CXX=g++
F90=gfortran
F77=gfortran

# cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
#         -DMPI_CXX_COMPILER=mpicxx \
#         -DLAPACK_DIR=$LAPACK \
#         -DSCALAPACK_DIR=$SCALAPACK \
#         -DELPA_DIR=$ELPA \
#         -DFFTW3_DIR=$FFTW3 \
#         -DCEREAL_INCLUDE_DIR=$CEREAL \
#         -DLibxc_DIR=$LIBXC \
#         -DENABLE_LCAO=ON \
#         -DENABLE_LIBXC=ON \
#         -DUSE_OPENMP=ON \
#         -DENABLE_ASAN=OFF \
#         -DUSE_ELPA=ON \
#         | tee configure.log

# if use deepks and deepmd
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DMPI_CXX_COMPILER=mpicxx \
        -DLAPACK_DIR=$LAPACK \
        -DSCALAPACK_DIR=$SCALAPACK \
        -DELPA_DIR=$ELPA \
        -DFFTW3_DIR=$FFTW3 \
        -DCEREAL_INCLUDE_DIR=$CEREAL \
        -DLibxc_DIR=$LIBXC \
        -DENABLE_LCAO=ON \
        -DENABLE_LIBXC=ON \
        -DUSE_OPENMP=ON \
        -DENABLE_ASAN=OFF \
        -DUSE_ELPA=ON \
        -DENABLE_DEEPKS=1 \
        -DTorch_DIR=$LIBTORCH \
        -Dlibnpy_INCLUDE_DIR=$LIBNPY \
	    -DDeePMD_DIR=$DEEPMD \
	    -DTensorFlow_DIR=$DEEPMD \
        | tee configure.log
# add mkl env for libtorch to link
module load mkl

cmake --build $BUILD_DIR -j `nproc` | tee build.log
cmake --install $BUILD_DIR | tee install.log
