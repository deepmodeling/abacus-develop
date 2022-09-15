#!/bin/bash
#
#PBS -l nodes=1:ppn=8
#PBS -l walltime=196:00:00
#PBS -j oe
#PBS -V
ulimit -s unlimited
cd $PBS_O_WORKDIR
NP=`cat $PBS_NODEFILE | wc -l`

cp INPUT-scf INPUT
cp KPT-scf KPT

mpirun -machinefile $PBS_NODEFILE -np $NP /public/udata/jingan/abacus_dft/abacus-NewGit/ABACUS.1.0.0/bin/ABACUS.mpi.1.0.0 > Log-scf.txt

cp INPUT-nscf INPUT
cp KPT-nscf KPT

mpirun -machinefile $PBS_NODEFILE -np $NP /public/udata/jingan/abacus_dft/abacus-NewGit/ABACUS.1.0.0/bin/ABACUS.mpi.1.0.0 > Log-nscf.txt

rm INPUT
rm KPT
