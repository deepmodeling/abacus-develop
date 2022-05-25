FROM ubuntu:latest
RUN apt update && apt install -y --no-install-recommends libopenblas-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc time sudo vim git # libomp-dev clang
RUN GIT_SSL_NO_VERIFY=true git clone https://gitee.com/deepmodeling/abacus-develop.git --depth 1 && cd abacus-develop && cmake -B build && cmake --build build -j8 && cmake --install build && cd .. && rm -rf abacus-develop
ENTRYPOINT abacus
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
