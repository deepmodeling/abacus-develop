#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 2 3 4; do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
    mpirun -np $i ./test_overlap_new_cd
    mpirun -np $i ./test_overlap_new
    mpirun -np $i ./test_ekinetic_new
    mpirun -np $i ./test_nonlocal_new
done
