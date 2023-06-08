#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp INPUT_1 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf1.output
mv OUT.autotest/running_scf.log OUT.autotest/running_scf1.log
cp INPUT_2 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf2.output
mv OUT.autotest/running_scf.log OUT.autotest/running_scf2.log

rm INPUT

if [[ ! -f scf1.output ]] ||
   [[ ! -f scf2.output ]] || 
   [[ ! -f OUT.autotest/running_scf1.log ]] ||
   [[ ! -f OUT.autotest/running_scf2.log ]] ||
   [[ ! ( "$(tail -1 OUT.autotest/running_scf1.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.autotest/running_scf2.log)" == " Total  Time  :"* ) ]] 

then
	echo "job is failed!"
	exit 1
else
	echo "job is successed!"
	exit 0
fi