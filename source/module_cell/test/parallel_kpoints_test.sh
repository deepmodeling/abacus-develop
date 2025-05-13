#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for ((i=3;i<=4;i++));
do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
    mpirun -np $i ./cell_ParaKpoints
    if [ $? -ne 0 ]; then
        echo "TEST in parallel, nprocs=$i failed"
        exit 1
    fi  
done
