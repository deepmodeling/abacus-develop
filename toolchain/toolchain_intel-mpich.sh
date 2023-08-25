#!/bin/bash
# JamesMisaka in 2023-08-25
# install abacus by intel-toolchain 
# use mkl and mpich, which is the fastest in test
# libtorch and libnpy are for deepks support, which can be =no
# can support deepmd 

# module load mkl mpi icc compiler

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