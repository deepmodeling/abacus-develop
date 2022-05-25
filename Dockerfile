# To build this docker file, run `docker build -t abacus - < Dockerfile`.
FROM ubuntu:latest
RUN apt update && apt install -y --no-install-recommends libopenblas-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc time sudo vim git
# If you wish to use the LLVM compiler, replace 'g++' above with 'clang libomp-dev'.
RUN GIT_SSL_NO_VERIFY=true git clone https://github.com/deepmodeling/abacus-develop.git --depth 1 && cd abacus-develop && cmake -B build && cmake --build build -j`nproc` && cmake --install build && cd .. && rm -rf abacus-develop
# If you have trouble cloning repo, replace "github.com" with "gitee.com".
ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none
CMD mpirun --use-hwthread-cpus abacus

# To run ABACUS built by this image with all available threads, execute `docker run -v <host>:<wd> -w <wd/input> abacus`.
# Replace '<host>' with the path to all files(including pseudopotential files), '<wd>' with a path to working directory, and '<wd/input>' with the path to input folder(containing 'INPUT', 'STRU', etc.).
# e.g. after clone the repo to `$HOME` and pulled this image, execute `docker run -v ~/abacus-develop/tests/integrate:/workspace -w /workspace/101_PW_15_f_pseudopots abacus`.
# To run ABACUS with a given MPI process number, execute `docker run -v <host>:<wd> -w <wd/input> -it --entrypoint mpirun abacus -np <processes> abacus`. Note: the first "abacus" is the name of the image, the second "abacus" is the name of the executable file. Do not use '--cpus' flag of 'docker run' to specify the number of processes.

# To use this image as developing environment, execute `docker run -it --entrypoint /bin/bash abacus`.
# Please refer to https://docs.docker.com/engine/reference/commandline/run/ for more details.
