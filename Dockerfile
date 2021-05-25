FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

RUN export DEBIAN_FRONTEND=noninteractive && apt update && \
    apt install -y wget gnupg gcc g++ gfortran cmake pkg-config build-essential git && \
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    echo "deb https://apt.repos.intel.com/oneapi all main" | tee -a /etc/apt/sources.list && \
    apt update && \
    apt install -y intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic intel-oneapi-compiler-fortran intel-oneapi-mkl-devel intel-oneapi-mpi-devel

RUN cd /opt/intel/oneapi/ && \
    rm -r conda_channel debugger dev-utilities tbb compiler/2021.2.0/linux/lib mpi/2021.2.0/lib/release mpi/2021.2.0/lib/release_mt mpi/2021.2.0/lib/debug_mt && \
    cd compiler/2021.2.0/linux/bin && find . -maxdepth 1 -type f -delete && cd ../../../.. && \
    mkdir mkl/2021.2.0/lib/tmp && cd mkl/2021.2.0/lib/tmp && \
    mv ../intel64/* ./ && mv libmkl_avx2.so.1 libmkl_avx512.so.1 libmkl_core.so* libmkl_intel_lp64.so* libmkl_blacs_intelmpi_lp64.so libmkl_intel_thread.so* libmkl_scalapack_lp64.so* locale ../intel64/ && \
    cd ../../../.. && rm -r mkl/2021.2.0/lib/tmp
