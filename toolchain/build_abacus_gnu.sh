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
LAPACK=${OPENBLAS_ROOT}/lib
SCALAPACK=${SCALAPACK_ROOT}/lib
ELPA=${ELPA_ROOT}
FFTW3=${FFTW_ROOT}
CEREAL=${CEREAL_ROOT}/include
LIBXC=${LIBXC_ROOT}
RAPIDJSON=${RAPIDJSON_ROOT}
LIBRI=${LIBRI_ROOT}
LIBCOMM=${LIBCOMM_ROOT}
USE_CUDA=OFF  # set ON to enable gpu-abacus

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
