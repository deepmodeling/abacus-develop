#!/bin/bash
#SBATCH -J install
#SBATCH -N 1
#SBATCH -n 64
# JamesMisaka in 2023-08-31
# install abacus by gnu-toolchain
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
