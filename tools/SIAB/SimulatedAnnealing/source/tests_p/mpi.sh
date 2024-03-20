#!/bin/bash
#
#PBS -l nodes=2:ppn=4
#PBS -j oe
#PBS -V
#PBS -N MCSP

# go to work dir
cd $PBS_O_WORKDIR

# setup Nums of Processor
NP=`cat $PBS_NODEFILE|wc -l`
echo "Numbers of Processors:  $NP"
echo "---------------------------"

# running program
/share/data1/mpich-intel/bin/mpirun -np $NP -machinefile $PBS_NODEFILE ../040.MCSP_p.exe > Log.txt
