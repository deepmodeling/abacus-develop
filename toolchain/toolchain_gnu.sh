#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# JamesMisaka in 2023-09-16
# install abacus dependency by gnu-toolchain
# one can use mpich or openmpi.
# openmpi will be faster, but not compatible in some cases.
# libtorch and libnpy are for deepks support, which can be =no
# if you want to run EXX calculation, you should set --with-libri=install
# mpich (and intel toolchain) is recommended for EXX support
# gpu-lcao supporting modify: CUDA_PATH and --enable-cuda
# export CUDA_PATH=/usr/local/cuda

./install_abacus_toolchain.sh \
--with-gcc=system \
--with-intel=no \
--with-openblas=install \
--with-openmpi=install \
--with-cmake=install \
--with-scalapack=install \
--with-libxc=install \
--with-fftw=install \
--with-elpa=install \
--with-cereal=install \
--with-rapidjson=install \
--with-libtorch=no \
--with-libnpy=no \
--with-libri=no \
--with-libcomm=no \
--with-4th-openmpi=no \
| tee compile.log
# to use openmpi-version4: set --with-4th-openmpi=yes
# to enable gpu-lcao, add the following lines:
# --enable-cuda \
# --gpu-ver=75 \ 
# one should check your gpu compute capability number 
# and use it in --gpu-ver
