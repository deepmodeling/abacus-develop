#!/bin/bash -e
#SBATCH -J build_abacus_gnu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err

# Build ABACUS by gnu-toolchain

# load system env modules at first

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)

BUILD_DIR=build_abacus_gnu
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
LAPACK=$INSTALL_DIR/openblas-0.3.30/lib
SCALAPACK=$INSTALL_DIR/scalapack-2.2.2/lib
ELPA=$INSTALL_DIR/elpa-2025.06.001/cpu
# ELPA=$INSTALL_DIR/elpa-2025.06.001/nvidia # for elpa-gpu
FFTW3=$INSTALL_DIR/fftw-3.3.10
CEREAL=$INSTALL_DIR/cereal-master/include/cereal
LIBXC=$INSTALL_DIR/libxc-7.0.0
RAPIDJSON=$INSTALL_DIR/rapidjson-master/
LIBRI=$INSTALL_DIR/LibRI-master
LIBCOMM=$INSTALL_DIR/LibComm-master
USE_CUDA=OFF  # set ON to enable gpu-abacus
# NEP_DIR=$INSTALL_DIR/NEP_CPU-main
# LIBTORCH=$INSTALL_DIR/libtorch-2.1.2/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-1.0.1/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd #

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
        -DENABLE_RAPIDJSON=ON \
        -DRapidJSON_DIR=$RAPIDJSON \
        -DENABLE_LIBRI=ON \
        -DLIBRI_DIR=$LIBRI \
        -DLIBCOMM_DIR=$LIBCOMM \
        -DUSE_CUDA=$USE_CUDA \
#         -DCMAKE_CUDA_COMPILER=/path/to/cuda/bin/nvcc \
#         -DNEP_DIR=$NEP_DIR \
#         -DENABLE_MLALGO=1 \
#         -DTorch_DIR=$LIBTORCH \
#         -Dlibnpy_INCLUDE_DIR=$LIBNPY \
# 	      -DDeePMD_DIR=$DEEPMD \
#         -DENABLE_CUSOLVERMP=ON \
#         -D CAL_CUSOLVERMP_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/2x.xx/math_libs/1x.x/targets/x86_64-linux/lib

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
EOF
