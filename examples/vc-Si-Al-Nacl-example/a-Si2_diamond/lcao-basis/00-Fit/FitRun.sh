#!/bin/bash

stru_start=10.00
max=40
delta=0.01
count=1
#element="Si"

while [ $count -le $max ]
do

latA=`echo " $stru_start + $count*$delta" |bc`
stru=`echo " ($stru_start + $count*$delta)* 1.88972687777" |bc`
echo " stru_num=stru-$latA"
#relax_nmax=1 

file="out_data_$latA"
test -d $file || mkdir $file

#mkdir $dis
cd $file 
cat > STRU << EOF
ATOMIC_SPECIES
Si 1.000 Si.pz-vbc.UPF 

NUMERICAL_ORBITAL
Si_50Ry_8.0au_dzp

LATTICE_CONSTANT
$latA

LATTICE_VECTORS
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5

ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.

Si      // Element type
0.0     // magnetism
2       // number of atoms
0.00 0.00 0.00 1 1 1
0.25 0.25 0.25 1 1 1

EOF
cp ../INPUT ./
cp ../KPT ./
cp ../Si.pz-vbc.UPF ./
cp ../Si_50Ry_8.0au_dzp ./
cp ../Submit_Einstein_intel.sh ./
#../RunTest 12
#cp ../node_openmpi .
#mpiexec -np 8 -machinefile node_openmpi ../../../bin/MESIA-ALPHA.fp_mpi.x 
#../../../bin/MESIA-ALPHA.fp_mpi.x
#mpiexec -n 6 -machinefile  ../../../bin/MESIA-ALPHA.fp_mpi.x > Log.txt
#/home/liuxiaohui/svn_test/bin/MESIA-ALPHA.fp_mpi.x > Log.txt $
qsub Submit_Einstein_intel.sh
cd ../ 

let count++
done
