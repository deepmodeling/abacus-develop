#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err
# install ABACUS with libxc and deepks
# JamesMisaka in 2023.08.31

# Build ABACUS by gnu-toolchain

# module load openmpi

TOOL=$(pwd)
ABACUS_DIR=..
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)

BUILD_DIR=build_abacus_gnu
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
LAPACK=$INSTALL_DIR/openblas-0.3.24/lib
SCALAPACK=$INSTALL_DIR/scalapalack-2.2.1/lib
ELPA=$INSTALL_DIR/elpa-2023.05.001/cpu
FFTW3=$INSTALL_DIR/fftw-3.3.10
CEREAL=$INSTALL_DIR/cereal-1.3.2/include/cereal
LIBXC=$INSTALL_DIR/libxc-6.2.2
# LIBRI=$INSTALL_DIR/LibRI-0.1.0
# LIBCOMM=$INSTALL_DIR/LibComm-0.1.0
# LIBTORCH=$INSTALL_DIR/libtorch-2.0.1/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-0.1.0/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=g++ \
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
        -DUSE_ELPA=ON \
#         -DENABLE_DEEPKS=1 \
#         -DTorch_DIR=$LIBTORCH \
#         -Dlibnpy_INCLUDE_DIR=$LIBNPY \
#         -DENABLE_LIBRI=ON \
#         -DLIBRI_DIR=$LIBRI \
#         -DLIBCOMM_DIR=$LIBCOMM \
# 	      -DDeePMD_DIR=$DEEPMD \
# 	      -DTensorFlow_DIR=$DEEPMD \

# # add mkl env for libtorch to link
# if one want to install libtorch, mkl should be load in build process
# for -lmkl when load libtorch
# module load mkl

# if one want's to include deepmd, your gcc version should be >= 11.3.0

cmake --build $BUILD_DIR -j `nproc` 
cmake --install $BUILD_DIR 2>/dev/null

# generate abacus_env.sh
cat << EOF > "${TOOL}/abacus_env.sh"
source $INSTALL_DIR/setup
export PATH="${PREFIX}/bin":${PATH}
EOF
