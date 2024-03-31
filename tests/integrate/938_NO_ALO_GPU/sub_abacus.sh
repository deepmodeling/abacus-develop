#!/bin/sh
#An example for MPI job.
#SBATCH -J abacus

#SBATCH -p normal
##SBATCH -p debug

#SBATCH -N 1 -n 32
##SBATCH --time=10-00

#SBATCH -o job-%j.log
#SBATCH -e job-%j.err

##SBATCH -n 512
##SBATCH --ntasks-per-node=8
##SBATCH --exclusive

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

module purge
module load compiler/devtoolset/9.3.1 compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export OMP_NUM_THREADS=1

MPIRUN=mpirun
#$MPIRUN /public/home/liuxh/software/code/abacus-develop-3.4.0/build-cpu/abacus
$MPIRUN /public/home/liuxh/software/code/abacus-develop-3.5.1/build-cpu/abacus
