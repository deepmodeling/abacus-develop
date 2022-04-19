#!/bin/bash

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 6 3 2;do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST LCAO DIAGO in parallel, nprocs=$i"
    mpirun -np $i ./hsolver_diago
    break    
done


