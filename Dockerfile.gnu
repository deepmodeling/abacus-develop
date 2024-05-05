# To build this Dockerfile, run `docker build -t abacus - < Dockerfile.gnu`.
# Build without cloning the repo by `docker build https://github.com/deepmodeling/abacus-develop.git#develop -f Dockerfile.intel`,
# Alternatively, pull the image with `docker pull ghcr.io/deepmodeling/abacus-intel:latest`.
# Also available at `docker pull registry.dp.tech/deepmodeling/abacus-intel:latest`.
# Available image names: abacus-gnu, abacus-intel, abacus-cuda

# Docker images are aimed for evaluating ABACUS.
# For production use, please compile ABACUS from source code and run in bare-metal for a better performace.

FROM ubuntu:22.04
RUN apt update && apt install -y --no-install-recommends \
    libopenblas-openmp-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev \
    libxc-dev libgtest-dev libgmock-dev libbenchmark-dev python3-numpy \
    bc cmake git g++ make bc time sudo unzip vim wget gfortran
    # If you wish to use the LLVM compiler, replace 'g++' above with 'clang libomp-dev'.

ENV GIT_SSL_NO_VERIFY=true TERM=xterm-256color \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none
    # The above environment variables are for using OpenMPI in Docker.

RUN git clone https://github.com/llohse/libnpy.git && \
    cp libnpy/include/npy.hpp /usr/local/include && \
    rm -r libnpy

RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip \
        --no-check-certificate --quiet -O libtorch.zip && \
    unzip -q libtorch.zip -d /opt && rm libtorch.zip

ENV CMAKE_PREFIX_PATH=/opt/libtorch/share/cmake

ADD https://api.github.com/repos/deepmodeling/abacus-develop/git/refs/heads/develop /dev/null
    # This will fetch the latest commit info, and store in docker building cache.
    # If there are newer commits, docker build will ignore the cache and build latest codes.

RUN git clone https://github.com/deepmodeling/abacus-develop.git --depth 1 && \
    cd abacus-develop && \
    cmake -B build -DENABLE_DEEPKS=ON -DENABLE_LIBXC=ON -DENABLE_LIBRI=ON -DENABLE_RAPIDJSON=ON && \
    cmake --build build -j`nproc` && \
    cmake --install build && \
    rm -rf build
    #&& rm -rf abacus-develop
