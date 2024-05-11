#!/bin/bash

CMAKE_INSTALL_PREFIX=/home/zaotian/Abacus/bin
CMAKE_CXX_COMPILER=mpiicpc
MPI_CXX_COMPILER=mpiicpc
ELPA_DIR=/home/software_app/elpa_lib_2021.05.002
ELPA_INCLUDE_DIRS=/home/software_app/elpa_lib_2021.05.002/include/elpa-2021.05.002
CEREAL_INCLUDE_DIR=/home/software_app/cereal-1.3.2

cmake -B build \
-DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
-DMPI_CXX_COMPILER=$MPI_CXX_COMPILER \
-DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
-DCEREAL_INCLUDE_DIR=$CEREAL_INCLUDE_DIR \
-DDEBUG_INFO=0 \
-DUSE_ELPA=0 \

cmake --build build/ -j 20
cmake --install build/
