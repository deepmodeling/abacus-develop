#!/bin/bash

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}'`
echo "nprocs in this machine is $np"

if [[ 2 -gt $np ]]; then
    echo "skip LOBPCG band-parallel UT: fewer than 2 cpu cores"
    exit 0
fi

echo "TEST DIAGO LOBPCG in band parallel, nprocs=2"
ABACUS_LOBPCG_TEST_BNDPAR=1 OMP_NUM_THREADS=1 mpirun -np 2 ./MODULE_HSOLVER_lobpcg \
    --gtest_filter=DiagoLobpcgTest.GeneralizedBandParallelRankCompressedSubspace:DiagoLobpcgTest.BandParallelReusesProjectedSearchDirectionProducts
e1=$?
if [[ e1 -ne 0 ]]; then
    echo -e "\e[1;33m [  FAILED  ] \e[0m execute LOBPCG band-parallel UT with 2 cores error."
    exit 1
fi
