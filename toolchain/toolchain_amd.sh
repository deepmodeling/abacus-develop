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

./install_abacus_toolchain.sh \
--with-amd=system \
--math-mode=aocl \
--with-intel=no \
--with-gcc=no \
--with-openmpi=install \
--with-cmake=install \
--with-scalapack=system \
--with-libxc=install \
--with-fftw=system \
--with-elpa=install \
--with-cereal=install \
--with-rapidjson=install \
--with-libtorch=no \
--with-libnpy=no \
--with-libri=no \
--with-libcomm=no \
--with-4th-openmpi=no \
--with-flang=no \
| tee compile.log
# to use openmpi-version4: set --with-4th-openmpi=yes
# flang is not recommended to use in this stage