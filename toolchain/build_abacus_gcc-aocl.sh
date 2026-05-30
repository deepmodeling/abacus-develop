#!/bin/bash -e
#SBATCH -J build_abacus_aocl
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o install.log
#SBATCH -e install.err

# Build ABACUS by gcc-aocl toolchain

# load system env modules at first
# module load openmpi aocc aocl

ABACUS_DIR=..
TOOL=$(pwd)
INSTALL_DIR=$TOOL/install
source $INSTALL_DIR/setup
cd $ABACUS_DIR
ABACUS_DIR=$(pwd)
#AOCLhome=/opt/aocl-linux-aocc-5.0.0/5.0.0/aocl/  # user should specify this parameter to the aocl installation path

BUILD_DIR=build_abacus_gcc_aocl
rm -rf $BUILD_DIR

PREFIX=$ABACUS_DIR
ELPA=${ELPA_ROOT}
CEREAL=${CEREAL_ROOT}/include
LIBXC=${LIBXC_ROOT}
RAPIDJSON=${RAPIDJSON_ROOT}
LAPACK=$AOCLhome/lib
SCALAPACK=$AOCLhome/lib
FFTW3=$AOCLhome
LIBRI=${LIBRI_ROOT}
LIBCOMM=${LIBCOMM_ROOT}
USE_CUDA=OFF  # set ON to enable gpu-abacus
# NEP_DIR=$INSTALL_DIR/NEP_CPU-main
# LIBTORCH=$INSTALL_DIR/libtorch-2.1.2/share/cmake/Torch
# LIBNPY=$INSTALL_DIR/libnpy-1.0.1/include
# DEEPMD=$HOME/apps/anaconda3/envs/deepmd 

NUM_JOBS="$(nproc)"
while [[ $# -gt 0 ]]; do
  case $1 in
    -j)
      if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then
        NUM_JOBS="${2}"
        shift 2
      else
        echo "ERROR: -j requires a number argument"
        exit 1
      fi
      ;;
    -j[0-9]*)
      NUM_JOBS="${1#-j}"
      shift
      ;;
    *)
      echo "ERROR: Unsupported argument: $1" >&2
      echo "Usage: $0 [-j N|-jN]" >&2
      exit 1
      ;;
  esac
done

cmake -B $BUILD_DIR -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_CXX_COMPILER=g++ \
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

cmake --build $BUILD_DIR --target install -j "${NUM_JOBS}"

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
