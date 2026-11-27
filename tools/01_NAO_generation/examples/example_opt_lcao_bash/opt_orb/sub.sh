#!/bin/sh
#BSUB -q renxg
#BSUB -oo job.log
#BSUB -n 24
#BSUB -m node151
export OMP_NUM_THREADS=8
EXEC=/home/nic/wszhang/eclipse_project/abacus-NewGit/ABACUS.1.0.0/tools/tools_exx/opt_orb_pytorch/main.py
python3 $EXEC > running_stdout.log
