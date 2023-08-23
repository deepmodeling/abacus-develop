#!/bin/bash
# JamesMisaka in 2023-08-22
# install abacus by intel-toolchain
# use mkl and intelmpi
# but mpich and openmpi can also be tried
# libtorch and libnpy are for deepks support, which can be =no

# module load mkl mpi icc compiler

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