FROM ubuntu:22.04
RUN apt update && apt install -y --no-install-recommends \
    libopenblas-openmp-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev \
    libxc-dev libgtest-dev libgmock-dev python3-numpy \
    bc cmake git g++ make bc time sudo unzip vim wget gfortran libmpich-dev mpich

ENV GIT_SSL_NO_VERIFY=true TERM=xterm-256color \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none

RUN git clone https://github.com/llohse/libnpy.git && \
    cp libnpy/include/npy.hpp /usr/local/include && \
    rm -r libnpy

RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip \
        --no-check-certificate --quiet -O libtorch.zip && \
    unzip -q libtorch.zip -d /opt  && rm libtorch.zip

ENV CMAKE_PREFIX_PATH=/opt/libtorch/share/cmake

ADD https://api.github.com/repos/deepmodeling/abacus-develop/git/refs/heads/develop /dev/null

ENV CMAKE_Fortran_COMPILER=/usr/bin/mpifort

RUN git clone https://github.com/deepmodeling/abacus-develop.git --depth 1 && \
    cd abacus-develop && \
    cmake -B build -DENABLE_DEEPKS=ON -DENABLE_LIBXC=ON -DENABLE_LIBRI=ON && \
    cmake --build build -j`nproc` && \
    cmake --install build && \
    rm -rf build && \
    cd .. 
    #&& rm -rf abacus-develop
