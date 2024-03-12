#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

# use OMP but not MPI in HSE
allprocess=$((${ABACUS_NPROCS} * ${ABACUS_THREADS}))
OMP_NUM_THREADS=${allprocess} mpirun -np 1 ${ABACUS_PATH} | tee scf.output

if [[ ! -f scf.output ]] || 
   [[ ! -f OUT.ABACUS/running_scf.log ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf.log)" == " Total  Time  :"* ) ]] 
then
	echo "job is failed!"
	exit 1
else
	echo "job is successed!"
	exit 0
fi