#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 2 3 4; do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
    mpirun -np $i ./operator_overlap_cd_test
    mpirun -np $i ./operator_overlap_test
    mpirun -np $i ./operator_ekinetic_test
    mpirun -np $i ./operator_nonlocal_test
    mpirun -np $i ./operator_T_NL_cd_test
done
