#!/bin/sh
#/share/data2/openmpi-gnu/bin/mpirun -np $1 -machinefile node_openmpi ../042.MCSP_p.exe
/opt/openmpi-intel9/bin/mpirun -np $1 -machinefile node_openmpi ../042.MCSP_p.exe

