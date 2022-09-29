#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp STRU_0 STRU
sed -i "/ntype/ c\ntype		2" INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H2O_scf.output
cp STRU_1 STRU
sed -i "/ntype/ c\ntype		2" INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee O_scf.output
cp STRU_2 STRU
sed -i "/ntype/ c\ntype		3" INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H1_scf.output
cp STRU_3 STRU
sed -i "/ntype/ c\ntype		3" INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H2_scf.output

rm STRU

E_H2O=`grep GE H2O_scf.output | tail -n 1 | awk '{printf $2}'`
E_O=`grep GE O_scf.output | tail -n 1 | awk '{printf $2}'`
E_H1=`grep GE H1_scf.output | tail -n 1 | awk '{printf $2}'`
E_H2=`grep GE H2_scf.output | tail -n 1 | awk '{printf $2}'`
awk 'BEGIN{print "'$E_H2O'" - "'$E_O'" - "'$E_H1'" - "'$E_H2'"}' > result.out
