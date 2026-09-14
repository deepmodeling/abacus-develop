#!/usr/bin/env bash

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export LD_PRELOAD=${LD_PRELOAD:-}
export CPATH=${CPATH:-}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

source /etc/profile.d/lmod.sh
module purge
module load cmake/3.31.6
module load abacus/develop-git-079fd0c-260724-sm70-auto

# CMake consumes the search paths exported by the loaded modules.
export CMAKE_LIBRARY_PATH=${LIBRARY_PATH:-}
export CMAKE_INCLUDE_PATH=${CPATH:-}
