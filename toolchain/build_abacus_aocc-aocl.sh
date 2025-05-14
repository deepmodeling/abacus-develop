#!/bin/bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err
# JamesMisaka in 2025.03.09

# Build ABACUS by amd-openmpi toolchain

# module load openmpi aocc aocl

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)
#AOCLhome=/opt/aocl-linux-aocc-5.0.0/5.0.0/aocl/  # user should specify this parameter

BUILD_DIR=build_abacus_aocl
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
ELPA=$INSTALL_DIR/elpa-2025.01.001/cpu
# ELPA=$INSTALL_DIR/elpa-2025.01.001/nvidia # for gpu-lcao
CEREAL=$INSTALL_DIR/cereal-master/include/cereal
LIBXC=$INSTALL_DIR/libxc-7.0.0
RAPIDJSON=$INSTALL_DIR/rapidjson-master/
LAPACK=$AOCLhome/lib
SCALAPACK=$AOCLhome/lib
FFTW3=$AOCLhome
# LIBRI=$INSTALL_DIR/LibRI-0.2.1.0
# LIBCOMM=$INSTALL_DIR/LibComm-master
# LIBTORCH=$INSTALL_DIR/libtorch-2.1.2/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-1.0.1/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd # v3.0 might have problem

# if clang++ have problem, switch back to g++

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DMPI_CXX_COMPILER=mpicxx \
        -DLAPACK_DIR=$LAPACK \
        -DSCALAPACK_DIR=$SCALAPACK \
        -DFFTW3_DIR=$FFTW3 \
        -DELPA_DIR=$ELPA \
        -DCEREAL_INCLUDE_DIR=$CEREAL \
        -DLibxc_DIR=$LIBXC \
        -DENABLE_LCAO=ON \
        -DENABLE_LIBXC=ON \
        -DUSE_OPENMP=ON \
        -DUSE_ELPA=ON \
        -DENABLE_RAPIDJSON=ON \
        -DRapidJSON_DIR=$RAPIDJSON \
#         -DENABLE_DEEPKS=1 \
#         -DTorch_DIR=$LIBTORCH \
#         -Dlibnpy_INCLUDE_DIR=$LIBNPY \
#         -DENABLE_LIBRI=ON \
#         -DLIBRI_DIR=$LIBRI \
#         -DLIBCOMM_DIR=$LIBCOMM \
# 	      -DDeePMD_DIR=$DEEPMD \
#         -DUSE_CUDA=ON \
#         -DENABLE_CUSOLVERMP=ON \
#         -D CAL_CUSOLVERMP_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/2x.xx/math_libs/1x.x/targets/x86_64-linux/lib

# if one want's to include deepmd, your system gcc version should be >= 11.3.0 for glibc requirements

cmake --build $BUILD_DIR -j `nproc` 
cmake --install $BUILD_DIR 2>/dev/null

# generate abacus_env.sh
cat << EOF > "${TOOL}/abacus_env.sh"
#!/bin/bash
source $INSTALL_DIR/setup
export PATH="${PREFIX}/bin":\${PATH}
EOF

# generate information
cat << EOF
========================== usage =========================
Done!
To use the installed ABACUS version
You need to source ${TOOL}/abacus_env.sh first !
"""
EOF