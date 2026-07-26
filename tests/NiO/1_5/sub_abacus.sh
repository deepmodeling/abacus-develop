#!/bin/bash
#SBATCH --job-name=NiO-mix-0.1
#SBATCH --partition=mars-gpu          
#SBATCH --nodes=1                     
#SBATCH --ntasks-per-node 2
#SBATCH --gres=gpu:2               
#SBATCH --output=abacus-spin.log           
#SBATCH --error=abacus-spin.err            
module load singularity
module load hpcc/ompi
module load hpcc/ucx
module list
srun --mpi=pmi2  singularity run --env OMP_NUM_THREADS=4  --env HPCC_CPU_THREAD_POLICY=2 --env HPCC_DIRECT_DISPATCH=0 -B /gpfs_ssd/metax:/gpfs_ssd/metax  --pwd ${PWD}  /gpfs_ssd/metax/software_sif/abacus-spin-u-1016.sif abacus_pw