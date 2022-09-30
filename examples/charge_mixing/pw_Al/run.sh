#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp INPUT_1 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf1.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf1.log
cp INPUT_2 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf2.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf2.log
cp INPUT_3 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf3.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf3.log
cp INPUT_4 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf4.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf4.log

rm INPUT
