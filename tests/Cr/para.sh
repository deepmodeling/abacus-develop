#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
source /public3/home/a6s000983/abacus-install/abacus-develop-mc/intel.env
export OMP_NUM_THREADS=64
export LD_LIBRARY_PATH=/public3/home/a6s000983/abacus-install/abacus-develop-mc/lib/:$LD_LIBRARY_PATH 
export LIBRARY_PATH=/public3/home/a6s000983/abacus-install/abacus-develop-mc/lib/:$LIBRARY_PATH
mpirun -n 1 /public3/home/a6s000983/abacus-install/abacus-develop-mc/intel/abacus
