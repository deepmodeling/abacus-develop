#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err
# install ABACUS with libxc and deepks
# JamesMisaka in 2023.08.22

# Build ABACUS by intel-toolchain

# module load mkl compiler mpi
# source path/to/vars.sh

TOOL=$(pwd)
ABACUS_DIR=..
BUILD_DIR=$ABACUS_DIR/build_abacus
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
rm -rf $BUILD_DIR
cd $ABACUS_DIR

PREFIX=$ABACUS_DIR
ELPA=$INSTALL_DIR/elpa-2021.11.002/cpu
CEREAL=$INSTALL_DIR/cereal-1.3.2/include/cereal
LIBXC=$INSTALL_DIR/libxc-6.2.2
LIBTORCH=$INSTALL_DIR/libtorch-2.0.1/share/cmake/Torch
LIBNPY=$INSTALL_DIR/libnpy-0.1.0/include
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
        -DENABLE_LIBRI=OFF \
        -DUSE_OPENMP=ON \
        -DENABLE_ASAN=OFF \
        -DUSE_ELPA=ON \
        -DENABLE_DEEPKS=1 \
        -DTorch_DIR=$LIBTORCH \
        -Dlibnpy_INCLUDE_DIR=$LIBNPY \
#         -DENABLE_LIBRI=ON \
#         -DLIBRI_DIR=$LIBRI \
#         -DLIBCOMM_DIR=$LIBCOMM \
# 	      -DDeePMD_DIR=$DEEPMD \
# 	      -DTensorFlow_DIR=$DEEPMD \

cmake --build $BUILD_DIR -j `nproc` 
cmake --install $BUILD_DIR 

# generate abacus_env.sh
cat << EOF > "${TOOL}/abacus_env.sh"
source $INSTALL_DIR/setup
export PATH="${PREFIX}":${PATH}
EOF
