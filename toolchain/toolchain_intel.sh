#!/bin/bash
# JamesMisaka in 2023-08-22
# install abacus by intel-toolchain
# use mkl and intelmpi
# but mpich and openmpi can also be tried

./install_abacus_toolchain.sh \
--with-intel=system --math-mode=mkl \
--with-gcc=no --with-intelmpi=system \
--with-cmake=install \
--with-scalapack=install \
--with-libxc=install \
--with-fftw=install \
--with-elpa=install \
--with-libtorch=no \
