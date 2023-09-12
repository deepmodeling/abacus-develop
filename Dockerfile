# To build this Dockerfile, run `docker build -t abacus - < Dockerfile`.
# Build without cloning the repo by `docker build https://github.com/deepmodeling/abacus-develop.git#develop`,
#   and optionally choose the Dockerfile in use by appending e.g. `-f Dockerfile.gnu`.
# Alternatively, pull the image with `docker pull ghcr.io/deepmodeling/abacus:latest`.
# Also available at `docker pull registry.dp.tech/deepmodeling/abacus:latest`.

# Docker images are aimed for evaluating ABACUS.
# For production use, please compile ABACUS from source code and run in bare-metal for a better performace.

FROM ubuntu:22.04
RUN apt update && apt install -y --no-install-recommends \
    libopenblas-openmp-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev \
    g++ make cmake bc time sudo vim git
# If you wish to use the LLVM compiler, replace 'g++' above with 'clang libomp-dev'.

ENV GIT_SSL_NO_VERIFY=true TERM=xterm-256color \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none

# This will fetch the latest commit info, and store in docker building cache.
# If there are newer commits, docker build will ignore the cache and build latest codes.
ADD https://api.github.com/repos/deepmodeling/abacus-develop/git/refs/heads/develop /dev/null

RUN git clone https://github.com/deepmodeling/abacus-develop.git --depth 1 && \
    cd abacus-develop && \
    cmake -B build && \
    cmake --build build -j`nproc` && \
    cmake --install build && \
    cd .. 
    #&& rm -rf abacus-develop
# If you have trouble cloning repo, replace "github.com" with "gitee.com".
CMD mpirun --use-hwthread-cpus abacus

# To run ABACUS built by this image with all available threads, execute `docker run -v <host>:<wd> -w <wd/input> abacus:latest`.
# Replace '<host>' with the path to all files(including pseudopotential files), '<wd>' with a path to working directory, and '<wd/input>' with the path to input folder(containing 'INPUT', 'STRU', etc.).
# e.g. after cloning the repo to `$HOME` and pulling image, execute `docker run -v ~/abacus-develop/tests/integrate:/workspace -w /workspace/101_PW_15_f_pseudopots abacus:latest`.
# To run ABACUS with a given MPI process number, execute `docker run -v <host>:<wd> -w <wd/input> -it --entrypoint mpirun abacus:latest -np <processes> abacus`.
# Note: It would be better using all available CPUs. Docker uses CFS to share the CPU resources, which will result in bad CPU affinity.

# To use this image as developing environment, execute `docker run -it --entrypoint /bin/bash abacus`.
# Please refer to https://docs.docker.com/engine/reference/commandline/run/ for more details.
