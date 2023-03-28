#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 4;do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
<<<<<<< HEAD
<<<<<<< HEAD:source/module_cell/test/bcast_atom_spec_test.sh
    mpirun -np $i ./cell_atom_spec
=======
    mpirun -np $i ./base_ParaReduce
>>>>>>> 597d101b5e2f0979645e60b803172ecac0895b52:source/module_base/test_parallel/parallel_reduce_test.sh
=======
    mpirun -np $i ./cell_atom_spec
>>>>>>> 597d101b5e2f0979645e60b803172ecac0895b52
    break    
done
