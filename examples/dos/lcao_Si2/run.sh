#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

cp INPUT1 INPUT
cp KPT1 KPT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee scf.output

cp INPUT2 INPUT
cp KPT2 KPT
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee nscf.output

rm INPUT KPT

if [[ ! -f scf.output ]] || 
   [[ ! -f nscf.output ]] ||
   [[ ! -f OUT.ABACUS/running_scf.log ]] ||
   [[ ! -f OUT.ABACUS/running_nscf.log ]] ||
   [[ ! -f OUT.ABACUS/DOS1 ]] ||
   [[ ! -f OUT.ABACUS/DOS1_smearing.dat ]] ||
   [[ ! -f OUT.ABACUS/SPIN1_CHG.cube ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_nscf.log)" == " Total  Time  :"* ) ]]
then
	echo "job is failed!"
	exit 1
else
	echo "job is successed!"
	exit 0
fi