#!/bin/bash
cmake -B build -DUSE_CUDA=ON -DDEBUG_INFO=ON -DUSE_CUSOLVER_LCAO=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j50
cmake --install build
cd examples/scf/lcao_Si2
cat /dev/null > ./result.txt
abacus >> result.txt