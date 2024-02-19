#!/bin/bash
cmake -B build -DUSE_CUDA=ON
cmake --build build -j50
cmake --install build
