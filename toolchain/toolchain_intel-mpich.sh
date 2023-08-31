#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 64

# JamesMisaka in 2023-08-25
# install abacus by intel-toolchain 
# use mkl , and mpich instead of intelmpi
# libtorch and libnpy are for deepks support, which can be =no
# can support deepmd 

# module load mkl compiler

./install_abacus_toolchain.sh \
--with-intel=system --math-mode=mkl \
--with-gcc=no --with-mpich=install \
--with-cmake=install \
--with-scalapack=no \
--with-libxc=install \
--with-fftw=no \
--with-elpa=install \
--with-cereal=install \
--with-libtorch=install \
--with-libnpy=install \
--with-intel-classic=yes \
| tee compile.log