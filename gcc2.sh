#!/bin/bash
> abacus.txt
cmake -B build -DCMAKE_CXX_COMPILER=/usr/bin/g++  \
               -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx \
               -DENABLE_LCAO=ON \
               -DENABLE_PAW=ON \
               -DENABLE_LIBXC=ON \
               -DENABLE_LIBRI=ON \
               -DENABLE_DEEPKS=ON \
               -DCMAKE_PREFIX_PATH=/home/ubuntu/anaconda3/envs/pytorch/lib/python3.10/site-packages/torch \
               -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.4/bin/nvcc \
               -DUSE_CUDA=ON 
cmake --build build -j50
cd /home/ubuntu/desktop/github/abacus/abacus-develop/tests/integrate/286_NO_KP_CR_HSE
OMP_NUM_THREADS=4 mpirun -n 2 /home/ubuntu/desktop/github/abacus/abacus-develop/build/abacus >> \
   /home/ubuntu/desktop/github/abacus/abacus-develop/abacus.txt
