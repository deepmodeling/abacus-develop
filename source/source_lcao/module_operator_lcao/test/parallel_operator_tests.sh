#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in 2 3 4; do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
    mpirun --allow-run-as-root -np $i ./MODULE_LCAO_operator_overlap_cd_test
    e1=$?
    mpirun --allow-run-as-root -np $i ./MODULE_LCAO_operator_overlap_test
    e2=$?
    mpirun --allow-run-as-root -np $i ./MODULE_LCAO_operator_ekinetic_test
    e3=$?
    mpirun --allow-run-as-root -np $i ./MODULE_LCAO_operator_nonlocal_test
    e4=$?
    mpirun --allow-run-as-root -np $i ./MODULE_LCAO_operator_T_NL_cd_test
    e5=$?
    if [[ $e1 -ne 0 || $e2 -ne 0 || $e3 -ne 0 || $e4 -ne 0 || $e5 -ne 0 ]]; then
        echo -e "\e[1;33m [  FAILED  ] \e[0m"\
			"execute UT with $i cores error."
        exit 1
    fi
done
