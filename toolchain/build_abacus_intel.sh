#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o build_abacus.log
#SBATCH -e build_abacus.err
# install ABACUS with libxc and deepks
# JamesMisaka in 2023.08.22

# Build ABACUS by intel-toolchain

#rm -rf ../build_abacus
#rm -rf ../build_abacus
# module load mkl compiler mpi
# source path/to/vars.sh

TOOL=$(pwd)
ABACUS_DIR=..
source ./install/setup
cd $ABACUS_DIR

PREFIX=.
BUILD_DIR=build_abacus
ELPA=$TOOL/install/elpa-2021.11.002/cpu
CEREAL=$TOOL/install/cereal-1.3.2/include/cereal
LIBXC=$TOOL/install/libxc-6.2.2
LIBTORCH=$TOOL/install/libtorch-2.0.1/share/cmake/Torch
LIBNPY=$TOOL/install/libnpy-0.1.0/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd

# if use deepks and deepmd
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=icpc \
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
        -DENABLE_DEEPKS=1 \
        -DTorch_DIR=$LIBTORCH \
        -Dlibnpy_INCLUDE_DIR=$LIBNPY \
        | tee configure.log
# 	      -DDeePMD_DIR=$DEEPMD \
# 	      -DTensorFlow_DIR=$DEEPMD \

cmake --build $BUILD_DIR -j `nproc` | tee build.log
cmake --install $BUILD_DIR | tee build.log
