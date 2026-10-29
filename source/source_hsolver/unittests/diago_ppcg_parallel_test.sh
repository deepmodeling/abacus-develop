#!/bin/bash

ppcg_binary=${1:-./MODULE_HSOLVER_ppcg_parallel}
np=$(nproc --all)
echo "nprocs in this machine is $np"

runs=0
for i in 2 3 6; do
    if [[ $i -gt $np ]]; then
        continue
    fi
    echo "TEST DIAGO PPCG in parallel, nprocs=$i"
    OMP_NUM_THREADS=1 mpirun -np $i "$ppcg_binary"
    e=$?
    if [[ $e -ne 0 ]]; then
        echo -e "\e[1;33m [  FAILED  ] \e[0m"\
            "execute UT with $i cores error."
        exit 1
    fi
    runs=$((runs + 1))
done

if [[ $np -ge 4 ]]; then
    echo "TEST DIAGO PPCG with two pools, nprocs=4"
    OMP_NUM_THREADS=1 mpirun -np 4 "$ppcg_binary" 2
    e=$?
    if [[ $e -ne 0 ]]; then
        echo -e "\e[1;33m [  FAILED  ] \e[0m"\
            "execute two-pool UT with 4 cores error."
        exit 1
    fi
    runs=$((runs + 1))
fi

if [[ $runs -eq 0 ]]; then
    echo "No 2/3/6-process case is supported; running the 1-process fallback."
    OMP_NUM_THREADS=1 mpirun -np 1 "$ppcg_binary"
    e=$?
    if [[ $e -ne 0 ]]; then
        echo -e "\e[1;33m [  FAILED  ] \e[0m"\
            "execute UT with 1 core error."
        exit 1
    fi
fi
