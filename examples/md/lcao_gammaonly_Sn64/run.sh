#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp INPUT_0 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md0.output
cp OUT.ABACUS/running_md.log running_md0.log
cp INPUT_1 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md1.output
cp OUT.ABACUS/running_md.log running_md1.log
cp INPUT_2 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md2.output
cp OUT.ABACUS/running_md.log running_md2.log
cp INPUT_3 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md3.output
cp OUT.ABACUS/running_md.log running_md3.log
cp INPUT_4 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md4.output
cp OUT.ABACUS/running_md.log running_md4.log

rm INPUT
