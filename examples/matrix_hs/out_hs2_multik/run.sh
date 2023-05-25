#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee output

if [[ ! -f output ]] || 
   [[ ! -f OUT.autotest/running_scf.log ]] ||
   [[ ! -f OUT.autotest/data-HR-sparse_SPIN0.csr ]] ||
   [[ ! -f OUT.autotest/data-SR-sparse_SPIN0.csr ]] ||
   [[ ! ( "$(tail -1 OUT.autotest/running_scf.log)" == " Total  Time  :"* ) ]] 
then
	echo "job is failed!"
	exit 1
else
	echo "job is successed!"
	exit 0
fi