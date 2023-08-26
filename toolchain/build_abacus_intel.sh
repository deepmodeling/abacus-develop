#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 64
# install ABACUS with libxc and deepks
# JamesMisaka in 2023.08.22

# Build ABACUS by intel-toolchain

#rm -rf ../build
# module load mkl mpi icc compiler
TOOL=$(pwd)
ABACUS_DIR=..
source ./install/setup
cd $ABACUS_DIR

PREFIX=./bin/abacus/
BUILD_DIR=build_abacus
ELPA=$TOOL/install/elpa-2021.11.002/cpu
CEREAL=$TOOL/install/cereal-1.3.2
LIBXC=$TOOL/install/libxc-6.2.2
# LIBTORCH=$TOOL/install/libtorch-2.0.1/share/cmake/Torch
# LIBNPY=$TOOL/install/libnpy-0.1.0/
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd

CC=icc
CXX=icpc
F90=ifort
F77=ifort

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DMPI_CXX_COMPILER=mpiicpc \
        -DMKLROOT=$MKLROOT \
        -DELPA_DIR=$ELPA \
        -DCEREAL_INCLUDE_DIR=$CEREAL \
        -DLibxc_DIR=$LIBXC \
        -DENABLE_LCAO=ON \
        -DENABLE_LIBXC=ON \
        -DUSE_OPENMP=ON \
        -DENABLE_ASAN=OFF \
        -DUSE_ELPA=ON \
        | tee configure.log

# if use deepks and deepmd
# cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
#         -DMPI_CXX_COMPILER=mpiicpc \
#         -DMKLROOT=$MKLROOT \
#         -DELPA_DIR=$ELPA \
#         -DCEREAL_INCLUDE_DIR=$CEREAL \
#         -DLibxc_DIR=$LIBXC \
#         -DENABLE_LCAO=ON \
#         -DENABLE_LIBXC=ON \
#         -DUSE_OPENMP=ON \
#         -DENABLE_ASAN=OFF \
#         -DUSE_ELPA=ON \
#         -DENABLE_DEEPKS=1 \
#         -DTorch_DIR=$LIBTORCH \
#         -Dlibnpy_INCLUDE_DIR=$LIBNPY \
# 	      -DDeePMD_DIR=$DEEPMD \
# 	      -DTensorFlow_DIR=$DEEPMD \
#         | tee configure.log

cmake --build $BUILD_DIR -j `nproc` | tee build.log
cmake --install $BUILD_DIR | tee install.log
