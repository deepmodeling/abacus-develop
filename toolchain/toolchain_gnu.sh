#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o compile.log
#SBATCH -e compile.err

# JamesMisaka in 2023-09-16
# install abacus dependency by gnu-toolchain
# one can use mpich or openmpi
# libtorch and libnpy are for deepks support, which can be =no

./install_abacus_toolchain.sh --with-openmpi=install \
--with-intel=no --with-gcc=system \
--with-cmake=install \
--with-scalapack=install \
--with-libxc=install \
--with-fftw=install \
--with-elpa=install \
--with-cereal=install \
--with-libtorch=no \
--with-libnpy=no \
--with-libri=no \
--with-libcomm=no \
| tee compile.log