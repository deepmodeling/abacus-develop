#!/bin/bash
#SBATCH -p cp6
#SBATCH -N 1
#SBATCH -J NaYbS2_LDA
#SBATCH -n 56

module add loginnode
source /fs2/software/intel/oneapi/2023.0/setvars.sh --force
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6
export PATH=$PATH:/fs2/home/chenmh/4_baotaoni/Softwares/abacus-develop/build
mpirun -n 9 abacus