#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

for i in $( seq 0 7 )
do
    cp INPUT_$i INPUT
    OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md$i.output
done

rm INPUT

outputs=(md0.output md1.output md2.output md3.output md4.output md5.output md6.output md7.output)
logs=(OUT.Si_nve/running_md.log OUT.Si_nhc_nvt/running_md.log OUT.Si_lgv/running_md.log OUT.Si_anderson/running_md.log OUT.Si_msst/running_md.log OUT.Si_berendsen/running_md.log OUT.Si_rescaling/running_md.log OUT.Si_rescale_v/running_md.log)

allpass=1

for i in "${outputs[@]}";do
        if [ ! -f "$i" ];then
                echo "No file: $i !"
                allpass=0
        fi
done

for i in "${logs[@]}";do
        if [ ! -f "$i" ];then
                echo "No file: $i !"
                allpass=0
        elif [[ ! ( "$(tail -1 $i)" == " Total  Time  :"* ) ]];then
                echo "File is not normal end: $i !"
                allpass=0
        fi
done

if [ $allpass == 0 ];then
        echo "job is failed!"
        exit 1
else
        echo "job is successed!"
        exit 0
fi
