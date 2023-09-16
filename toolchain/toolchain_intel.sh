#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# JamesMisaka in 2023-08-31
# install abacus by intel-toolchain
# use mkl and intelmpi
# but mpich and openmpi can also be tried
# libtorch and libnpy are for deepks support, which can be =no

# module load mkl mpi compiler

./install_abacus_toolchain.sh \
--with-intel=system --math-mode=mkl \
--with-gcc=no --with-intelmpi=system \
--with-cmake=install \
--with-scalapack=no \
--with-libxc=install \
--with-fftw=no \
--with-elpa=install \
--with-cereal=install \
--with-libtorch=install \
--with-libnpy=install \
--with-libri=install \
--with-libcomm=install \
--with-intel-classic=yes \
| tee compile.log