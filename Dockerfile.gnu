FROM debian:bullseye-slim

RUN apt-get update && apt-get install -y --no-install-recommends git g++ gfortran libssl-dev make cmake vim wget bc unzip \
    && apt-get install -y --no-install-recommends mpich libmpich-dev

ENV GIT_SSL_NO_VERIFY 1

RUN cd /tmp \
    && git clone https://github.com/USCiLab/cereal.git \
    && cp -r cereal/include /usr/local \
    && rm -rf cereal

RUN cd /tmp \
    && git clone https://github.com/xianyi/OpenBLAS.git --single-branch --depth=1 \
    && cd OpenBLAS && make USE_OPENMP=1 NO_AVX512=1 FC=gfortran -j8 && make PREFIX=/usr/local install \
    && cd /tmp && rm -rf OpenBLAS

RUN cd /tmp \
    && git clone https://github.com/darelbeida/scalapack.git -b v2.0.2-openblas --single-branch --depth=1 \
    && cd scalapack && make lib && cp libscalapack.a /usr/local/lib/ \
    && cd /tmp && rm -rf scalapack

RUN cd /tmp \
    && wget https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/2021.05.002/elpa-2021.05.002.tar.gz --no-check-certificate --quiet \
    && tar xzf elpa-2021.05.002.tar.gz && rm elpa-2021.05.002.tar.gz \
    && cd elpa-2021.05.002 && mkdir build && cd build \
    && ../configure CFLAGS="-O3 -march=native -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
    FCFLAGS="-O2 -mavx" --disable-avx512 \
    && make -j8 && make PREFIX=/usr/local install \
    && ln -s /usr/local/include/elpa-2021.05.002/elpa /usr/local/include/ \
    && cd /tmp && rm -rf elpa-2021.05.002

RUN cd /tmp \
    && wget http://www.fftw.org/fftw-3.3.9.tar.gz --no-check-certificate --quiet \
    && tar zxvf fftw-3.3.9.tar.gz \
    && cd fftw-3.3.9 \
    && ./configure --enable-mpi-fortran --enable-orterun-prefix-by-default FC=gfortran \
    && make -j8 && make PREFIX=/usr/local install \
    && cd /tmp && rm -rf fftw-3.3.9 && rm fftw-3.3.9.tar.gz

RUN cd /tmp \
    && wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.9.1%2Bcpu.zip --no-check-certificate --quiet \
    && unzip libtorch-cxx11-abi-shared-with-deps-1.9.1+cpu.zip \
    && cp -r libtorch/include /usr/local \
    && cp -r libtorch/lib /usr/local \
    && cp -r libtorch/share /usr/local \
    && rm -rf libtorch

RUN cd /tmp \
    && wget https://gitlab.com/libxc/libxc/-/archive/5.1.5/libxc-5.1.5.tar.gz --no-check-certificate --quiet \
    && tar xzf libxc-5.1.5.tar.gz \
    && cd libxc-5.1.5 \
    && mkdir build \
    && cmake -B build -DBUILD_TESTING=OFF \
    && cmake --build build \
    && cmake --install build \
    && cd /tmp \
    && rm -rf libxc-5.1.5 \
    && rm libxc-5.1.5.tar.gz

RUN cd /tmp \
    && git clone https://github.com/llohse/libnpy.git \
    && cp libnpy/include/npy.hpp /usr/local \
    && rm -rf libnpy

RUN cd /tmp \
    && git clone https://github.com/google/googletest.git \
    && cd googletest && cmake . && make install \
    && rm -rf googletest
