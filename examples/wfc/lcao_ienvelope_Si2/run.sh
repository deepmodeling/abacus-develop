#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp INPUT1 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf.output
cp INPUT2 INPUT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee get_wf.output

rm INPUT

if [[ ! -f scf.output ]] ||
   [[ ! -f get_wf.output ]] || 
   [[ ! -f OUT.ABACUS/running_scf.log ]] ||
   [[ ! -f OUT.ABACUS/running_get_wf.log ]] ||
   [[ ! -f OUT.ABACUS/WFC_NAO_K1.txt ]] ||
   [[ ! -f OUT.ABACUS/WFC_NAO_K36.txt ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_get_wf.log)" == " Total  Time  :"* ) ]] 
then
	echo "job failed!"
	exit 1
else
	echo "job succeeded!"
	exit 0
fi
