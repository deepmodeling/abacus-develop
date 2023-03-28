#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 4;do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
<<<<<<< HEAD
    mpirun -np $i ./Cell_ParaKpoints
=======
    mpirun -np $i ./cell_ParaKpoints
>>>>>>> 597d101b5e2f0979645e60b803172ecac0895b52
    break    
done
