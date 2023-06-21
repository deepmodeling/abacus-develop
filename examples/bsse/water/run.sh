#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

if [[ $ABACUS_NPROCS -gt 6 ]];then
	ABACUS_NPROCS=8
fi

cp STRU_0 STRU
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H2O_scf.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf_H2O.log
cp STRU_1 STRU
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee O_scf.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf_O.log
cp STRU_2 STRU
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H1_scf.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf_H1.log
cp STRU_3 STRU
OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee H2_scf.output
mv OUT.ABACUS/running_scf.log OUT.ABACUS/running_scf_H2.log

rm STRU

E_H2O=`grep GE H2O_scf.output | tail -n 1 | awk '{printf $2}'`
E_O=`grep GE O_scf.output | tail -n 1 | awk '{printf $2}'`
E_H1=`grep GE H1_scf.output | tail -n 1 | awk '{printf $2}'`
E_H2=`grep GE H2_scf.output | tail -n 1 | awk '{printf $2}'`
result=$(awk 'BEGIN{print "'$E_H2O'" - "'$E_O'" - "'$E_H1'" - "'$E_H2'"}')
echo $result > result.out
result_ref=$(cat result.ref)
difference=$(awk 'BEGIN{print "'$result'" - "'$result_ref'"}')
difference=$(awk -v a=$difference 'BEGIN{if(a<0) print -1*a; else print a}')

if [[ ! -f H2O_scf.output ]] || 
   [[ ! -f O_scf.output ]] ||
   [[ ! -f H1_scf.output ]] ||
   [[ ! -f H2_scf.output ]] ||
   [[ ! -f OUT.ABACUS/running_scf_H2O.log ]] ||
   [[ ! -f OUT.ABACUS/running_scf_O.log ]] ||
   [[ ! -f OUT.ABACUS/running_scf_H1.log ]] ||
   [[ ! -f OUT.ABACUS/running_scf_H2.log ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf_H2O.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf_O.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf_H1.log)" == " Total  Time  :"* ) ]] ||
   [[ ! ( "$(tail -1 OUT.ABACUS/running_scf_H2.log)" == " Total  Time  :"* ) ]] ||
   [[ $difference > 0.000001 ]] 
then
	echo "job is failed!"
	exit 1
else
	echo "job is successed!"
	exit 0
fi
