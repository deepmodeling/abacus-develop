#!/bin/bash
set -eu

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

np=$(awk '/cpu cores/ { print $NF; exit }' /proc/cpuinfo)
echo "nprocs in this machine is $np"

for ranks in 2 4; do
    if ((ranks > np)); then
        continue
    fi
    echo "TEST in parallel, nprocs=$ranks"
    if ! mpirun -np "$ranks" "${1:-./MODULE_LCAO_tddft_boundary_fix_test}"; then
        echo -e "\e[1;33m [  FAILED  ] \e[0m execute UT with $ranks ranks."
        exit 1
    fi
done
