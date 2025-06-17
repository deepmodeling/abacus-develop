#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

if [[ $ABACUS_NPROCS -gt 8 ]];then
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

E_H2O=$(grep FINAL_ETOT_IS OUT.ABACUS/running_scf_H2O.log | awk '{printf $2}')
E_O=$(grep FINAL_ETOT_IS OUT.ABACUS/running_scf_O.log | awk '{printf $2}')
E_H1=$(grep FINAL_ETOT_IS OUT.ABACUS/running_scf_H1.log | awk '{printf $2}')
E_H2=$(grep FINAL_ETOT_IS OUT.ABACUS/running_scf_H2.log | awk '{printf $2}')
result=$(echo "scale=12; ${E_H2O} - ${E_O} - ${E_H1} - ${E_H2}" | bc -l)
echo $result > result.out
echo "E_H2O: $E_H2O" >> result.out
echo "E_O: $E_O" >> result.out
echo "E_H1: $E_H1" >> result.out
echo "E_H2: $E_H2" >> result.out
result_ref=$(cat result.ref | head -1)
difference=$(echo "scale=12; $result - $result_ref" | bc -l)
abs_difference=$(echo "scale=12; if ($difference < 0) $difference * -1 else $difference" | bc)

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
   [ $(echo "$abs_difference < 0.00001" | bc) -ne 1 ]
then
	echo "job failed!"
	exit 1
else
	echo "job succeeded!"
	exit 0
fi
