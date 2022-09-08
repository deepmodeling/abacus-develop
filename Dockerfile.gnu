FROM ubuntu:22.04
RUN apt update && apt install -y --no-install-recommends \
    libopenblas-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev \
    libxc-dev libgtest-dev libgmock-dev python3-numpy \
    bc cmake git g++ make bc time sudo unzip vim wget

ENV GIT_SSL_NO_VERIFY=true TERM=xterm-256color \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none

RUN git clone https://github.com/llohse/libnpy.git && \
    cp libnpy/include/npy.hpp /usr/local/include && \
    rm -r libnpy

RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.9.1%2Bcpu.zip \
        --no-check-certificate --quiet -O libtorch.zip && \
    unzip -q libtorch.zip && rm libtorch.zip && \
    cd libtorch && cp -r . /usr/local && \
    cd .. && rm -r libtorch

RUN git clone https://github.com/deepmodeling/abacus-develop.git --depth 1 && \
    cd abacus-develop && \
    cmake -B build -DENABLE_DEEPKS=ON && \
    cmake --build build -j`nproc` && \
    cmake --install build && \
    cd .. && rm -rf abacus-develop
