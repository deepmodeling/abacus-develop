#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err
# build and install ABACUS with libxc, also can with deepks and deepmd
# JamesMisaka in 2023.08.31

# Build ABACUS by intel-toolchain with mpich

# module load mkl compiler
# source path/to/vars.sh

TOOL=$(pwd)
ABACUS_DIR=..
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)

BUILD_DIR=build_abacus
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
ELPA=$INSTALL_DIR/elpa-2023.05.001/cpu
CEREAL=$INSTALL_DIR/cereal-1.3.2/include/cereal
LIBXC=$INSTALL_DIR/libxc-6.2.2
# LIBTORCH=$INSTALL_DIR/libtorch-2.0.1/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-0.1.0/include
# LIBRI=$INSTALL_DIR/LibRI-0.1.0
# LIBCOMM=$INSTALL_DIR/LibComm-0.1.0
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=icpc \
        -DMPI_CXX_COMPILER=mpicxx \
        -DMKLROOT=$MKLROOT \
        -DELPA_DIR=$ELPA \
        -DCEREAL_INCLUDE_DIR=$CEREAL \
        -DLibxc_DIR=$LIBXC \
        -DENABLE_LCAO=ON \
        -DENABLE_LIBXC=ON \
        -DUSE_OPENMP=ON \
        -DENABLE_ASAN=OFF \
        -DUSE_ELPA=ON \
        # -DENABLE_DEEPKS=1 \
        # -DTorch_DIR=$LIBTORCH \
        # -Dlibnpy_INCLUDE_DIR=$LIBNPY \
        # -DENABLE_LIBRI=ON \
        # -DLIBRI_DIR=$LIBRI \
        # -DLIBCOMM_DIR=$LIBCOMM \
	    # -DDeePMD_DIR=$DEEPMD \
	    # -DTensorFlow_DIR=$DEEPMD \

# if one want's to include deepmd, your gcc version should be >= 11.3.0

cmake --build $BUILD_DIR -j `nproc` 
cmake --install $BUILD_DIR 

# generate abacus_env.sh
cat << EOF > "${TOOL}/abacus_env.sh"
source $INSTALL_DIR/setup
export PATH="${PREFIX}/bin":${PATH}
EOF