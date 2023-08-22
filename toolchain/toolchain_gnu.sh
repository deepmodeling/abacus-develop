#!/bin/bash
# JamesMisaka in 2023-08-22
# install abacus by gnu-toolchain
# one can use mpich or openmpi

./install_abacus_toolchain.sh --with-mpich=install \
--with-intel=no --with-gcc=system \
--with-cmake=install \
--with-scalapack=install \
--with-libxc=install \
--with-fftw=install \
--with-elpa=install \
--with-libtorch=no \
