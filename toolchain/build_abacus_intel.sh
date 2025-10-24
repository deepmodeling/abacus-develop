#!/bin/bash -e
#SBATCH -J build_abacus_intel
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err

# Build ABACUS by intel-toolchain

# load intel-oneapi env at first
# module load mkl compiler mpi
# source path/to/setvars.sh

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)

BUILD_DIR=build_abacus_intel
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
ELPA=$INSTALL_DIR/elpa-2024.05.001/cpu
# ELPA=$INSTALL_DIR/elpa-2024.05.001/nvidia # for gpu-lcao
CEREAL=$INSTALL_DIR/cereal-master/include
LIBXC=$INSTALL_DIR/libxc-7.0.0
RAPIDJSON=$INSTALL_DIR/rapidjson-master
LIBRI=$INSTALL_DIR/LibRI-master
LIBCOMM=$INSTALL_DIR/LibComm-master
USE_CUDA=OFF  # set ON to enable gpu-abacus
# NEP_DIR=$INSTALL_DIR/NEP_CPU-main
# LIBTORCH=$INSTALL_DIR/libtorch-2.1.2/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-1.0.1/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd # v3.0 might have problem

# Notice: if you are compiling with AMD-CPU or GPU-version ABACUS, then `icpc` and `mpiicpc` compilers are needed 
cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=icpx \
        -DMPI_CXX_COMPILER=mpiicpx \
        -DMKLROOT=$MKLROOT \
        -DELPA_DIR=$ELPA \
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
