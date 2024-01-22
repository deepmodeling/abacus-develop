# Easy Installation

This guide helps you install ABACUS with basic features. **For DeePKS, DeePMD and Libxc support, or building with `make`, please refer to [the advanced installation guide](../advanced/install.md)** after going through this page. We recommend building ABACUS with `cmake` to avoid dependency issues. We recommend compiling ABACUS(and possibly its requirements) from the source code using the latest compiler for the best performace. You can also deploy ABACUS **without building** by [Docker](#container-deployment) or [conda](#install-by-conda). Please note that ABACUS only supports Linux; for Windows users, please consider using [WSL](https://learn.microsoft.com/en-us/windows/wsl/) or docker.

## Prerequisites

To compile ABACUS, please make sure that the following prerequisites are present:

- [CMake](https://cmake.org/) >= 3.16 .
- C++ compiler, supporting C++11. You can use [Intel® C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/).

> GCC version 5 or later is always required. Intel compilers also use GCC headers and libraries[(ref)](https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compatibility-and-portability/gcc-compatibility-and-interoperability.html#gcc-compatibility-and-interoperability_GUID-52CB6FE0-83DA-4028-9EF4-0DFAF1652736).

- MPI library. The recommended versions are [Intel MPI](https://software.intel.com/enus/mpi-library), [MPICH](https://www.mpich.org/) or [Open MPI](https://www.open-mpi.org/).
- Fortran compiler if you are building `BLAS`, `LAPACK`, `ScaLAPACK`, and `ELPA` from source file. You can use [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or [GFortran](https://gcc.gnu.org/fortran/).
- [BLAS](http://www.netlib.org/blas/). You can use [OpenBLAS](https://www.openblas.net/).
- [LAPACK](http://www.netlib.org/lapack/).
- [FFTW3](http://www.fftw.org/).

These requirements support the calculation of plane-wave basis in ABACUS. For LCAO basis calculation, additional components are required:

- [ScaLAPACK](http://www.netlib.org/scalapack/).
- [CEREAL](https://uscilab.github.io/cereal/).
- [ELPA](https://elpa.mpcdf.mpg.de/) >= 2017 (optional).

## Install requirements

Some of these packages can be installed with popular package management system, such as `apt` and `yum`:

```bash
sudo apt update && sudo apt install -y libopenblas-openmp-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc git pkgconf
```

> Installing ELPA by apt only matches requirements on Ubuntu 22.04. For earlier linux distributions, you should build ELPA from source.

We recommend [Intel® oneAPI toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/commercial-base-hpc.html) (former Intel® Parallel Studio) as toolchain. The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#base-kit) contains Intel® oneAPI Math Kernel Library (aka `MKL`), including `BLAS`, `LAPACK`, `ScaLAPACK` and `FFTW3`. The [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#hpc-kit) contains Intel® MPI Library, and C++ compiler(including MPI compiler).
> Please note that building `elpa` with a different MPI library may cause conflict.
> Don't forget to [set environment variables](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-render-linux/top/configure-your-system.html) before you start! `cmake` will use Intel MKL if the environment variable `MKLROOT` is set.

Please refer to our [guide](https://github.com/deepmodeling/abacus-develop/wiki/Building-and-Running-ABACUS) on installing requirements.

## Install requirements by toolchain

We offer a set of [toolchain](https://github.com/deepmodeling/abacus-develop/tree/develop/toolchain)
scripts to compile and install all the requirements
automatically and suitable for machine characteristic in an online or offline way.
The toolchain can be downloaded with ABACUS repo, which is easily used and can
have a convenient installation under HPC environment in both `GNU` or `Intel-oneAPI` toolchain.
Sometimes, ABACUS by toolchain installation may have highly efficient performance.
A Tutorial for using this toolchain can be accessed in [bohrium-notebook](https://nb.bohrium.dp.tech/detail/5215742477)

> Notice: the toolchain is under development, please let me know if you encounter any problem in using this toolchain.


## Get ABACUS source code

Of course a copy of ABACUS source code is required, which can be obtained via one of the following choices:

- Clone the whole repo with git: `git clone https://github.com/deepmodeling/abacus-develop.git`
- Clone the minimum required part of repo: `git clone https://github.com/deepmodeling/abacus-develop.git --depth=1`
- Download the latest source code without git: `wget https://github.com/deepmodeling/abacus-develop/archive/refs/heads/develop.zip`
- Get the source code of a stable version [here](https://github.com/deepmodeling/abacus-develop/releases)
- If you have connection issues accessing GitHub, please try out our official [Gitee repo](https://gitee.com/deepmodeling/abacus-develop/): e.g. `git clone https://gitee.com/deepmodeling/abacus-develop.git`

### Update to latest release

Please check the [release page](https://github.com/deepmodeling/abacus-develop/releases) for the release note of a new version.

It is OK to download the new source code from beginning following the previous step.

To update your cloned git repo in-place:

```bash
git remote -v
# Check if the output contains the line below
# origin https://github.com/deepmodeling/abacus-develop.git (fetch)
# The remote name is marked as "upstream" if you clone the repo from your own fork.

# Replace "origin" with "upstream" or the remote name corresponding to deepmodeling/abacus-develop if necessary
git fetch origin
git checkout v3.2.0 # Replace the tag with the latest version
git describe --tags # Verify if the tag has been successfully checked out
```

Then proceed to the [Build and Install](#build-and-install) part. If you encountered errors, try remove the `build` directory first and reconfigure.

To use the codes under active development:

```bash
git checkout develop
git pull
```

## Configure

The basic command synopsis is:

```bash
cd abacus-develop
cmake -B build [-D <var>=<value>] ...
```

Here, 'build' is the path for building ABACUS; and '-D' is used for setting up some variables for CMake indicating optional components or requirement positions.

- `CMAKE_INSTALL_PREFIX`: the path of ABACUS binary to install; `/usr/local/bin/abacus` by default
- Compilers
  - `CMAKE_CXX_COMPILER`: C++ compiler; usually `g++`(GNU C++ compiler) or `icpx`(Intel C++ compiler). Can also set from environment variable `CXX`. It is OK to use MPI compiler here.
  - `MPI_CXX_COMPILER`: MPI wrapper for C++ compiler; usually `mpicxx` or `mpiicpc`(for Intel MPI).
- Requirements: Unless indicated, CMake will try to find under default paths.
  - `MKLROOT`: If environment variable `MKLROOT` exists, `cmake` will take MKL as a preference, i.e. not using `LAPACK`, `ScaLAPACK` and `FFTW`. To disable MKL, unset environment variable `MKLROOT`, or pass `-DMKLROOT=OFF` to `cmake`.
  - `LAPACK_DIR`: Path to OpenBLAS library `libopenblas.so`(including BLAS and LAPACK)
  - `SCALAPACK_DIR`: Path to ScaLAPACK library `libscalapack.so`
  - `ELPA_DIR`: Path to ELPA install directory; should be the folder containing 'include' and 'lib'.
  > Note: In ABACUS v3.5.1 or earlier, if you install ELPA from source , please add a symlink to avoid the additional include file folder with version name: `ln -s elpa/include/elpa-2021.05.002/elpa elpa/include/elpa` to help the build system find ELPA headers.

  - `FFTW3_DIR`: Path to FFTW3.
  - `CEREAL_INCLUDE_DIR`: Path to the parent folder of `cereal/cereal.hpp`. Will download from GitHub if absent.
  - `Libxc_DIR`: (Optional) Path to Libxc.
  > Note: In ABACUS v3.5.1 or earlier, Libxc built from source with Makefile is NOT supported; please compile Libxc with CMake instead.
  - `LIBRI_DIR`: (Optional) Path to LibRI.
  - `LIBCOMM_DIR`: (Optional) Path to LibComm.

- Components: The values of these variables should be 'ON', '1' or 'OFF', '0'. The default values are given below.
  - `ENABLE_LCAO=ON`: Enable LCAO calculation. If SCALAPACK, ELPA or CEREAL is absent and only require plane-wave calculations, the feature of calculating LCAO basis can be turned off.
  - `ENABLE_LIBXC=OFF`: [Enable Libxc](../advanced/install.md#add-libxc-support) to suppport variety of functionals. If `Libxc_DIR` is defined, `ENABLE_LIBXC` will set to 'ON'.
  - `ENABLE_LIBRI=OFF`: [Enable LibRI](../advanced/install.md#add-libri-support) to suppport variety of functionals. If `LIBRI_DIR` and `LIBCOMM_DIR` is defined, `ENABLE_LIBRI` will set to 'ON'.
  - `USE_OPENMP=ON`: Enable OpenMP support. Building ABACUS without OpenMP is not fully tested yet.
  - `BUILD_TESTING=OFF`: [Build unit tests](../advanced/install.md#build-unit-tests).
  - `ENABLE_MPI=ON`: Enable MPI parallel compilation. If set to `OFF`, a serial version of ABACUS with PW basis only will be compiled. Currently serial version of ABACUS with LCAO basis is not supported yet, so `ENABLE_LCAO` will be automatically set to `OFF`.
  - `ENABLE_COVERAGE=OFF`: Build ABACUS executable supporting [coverage analysis](../CONTRIBUTING.md#generating-code-coverage-report). This feature has a drastic impact on performance.
  - `ENABLE_ASAN=OFF`: Build with Address Sanitizer. This feature would help detecting memory problems.
  - `USE_ELPA=ON`: Use ELPA library in LCAO calculations. If this value is set to OFF, ABACUS can be compiled without ELPA library.

Here is an example:

```bash
CXX=mpiicpc cmake -B build -DCMAKE_INSTALL_PREFIX=~/abacus -DELPA_DIR=~/elpa-2016.05.004/build -DCEREAL_INCLUDE_DIR=~/cereal/include
```

## Build and Install

After configuring, build and install by:

```bash
cmake --build build -j`nproc`
cmake --install build
```

You can change the number after `-j` on your need: set to the number of CPU cores(`nproc`) to reduce compilation time.

## Run

If ABACUS is installed into a custom directory using `CMAKE_INSTALL_PREFIX`, please add it to your environment variable `PATH` to locate the correct executable.

```bash
export PATH=/my-install-dir/:$PATH
```

Please set OpenMP threads by setting environment variable:

```bash
export OMP_NUM_THREADS=1
```

Enter a directory containing a `INPUT` file. Please make sure structure, pseudo potential, or orbital files indicated by `INPUT` is at the correct location.

```bash
cd abacus-develop/examples/force/pw_Si2
```

Use 4 MPI processes to run, for example:

```bash
mpirun -n 4 abacus
```

> The total thread count(i.e. OpenMP per-process thread count * MPI process count) should not exceed the number of cores in your machine.

Please refer to [hands-on guide](./hands_on.md) for more instructions.

> Note: Some Intel CPU has a feature named Hyper-Threading(HT). This feature enables one physical core switch fastly between two logical threads. It would benefits from I/O bound tasks: when a thread is blocked by I/O, the CPU core can work on another thread. However, it helps little on CPU bound tasks, like ABACUS and many other scientific computing softwares. We recommend using the physical CPU core number.
> To determine if HT is turned on, execute `lscpu | grep 'per core'` and see if 'Thread(s) per core' is 2.

## Container Deployment

> Please note that containers target at developing and testing, but not massively parallel computing for production. Docker has a bad support to MPI, which may cause performance degradation.

We've built a ready-for-use version of ABACUS with docker [here](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus). For a quick start: pull the image, prepare the data, run container. Instructions on using the image can be accessed in [Dockerfile](../../Dockerfile). A mirror is available by `docker pull registry.dp.tech/deepmodeling/abacus`.

We also offer a pre-built docker image containing all the requirements for development. Please refer to our [Package Page](https://github.com/orgs/deepmodeling/packages?repo_name=abacus-develop).

The project is ready for VS Code development container. Please refer to [Developing inside a Container](https://code.visualstudio.com/docs/remote/containers#_quick-start-try-a-development-container). Choose `Open a Remote Window -> Clone a Repository in Container Volume` in VS Code command palette, and put the [git address](https://github.com/deepmodeling/abacus-develop.git) of `ABACUS` when prompted.

For online development environment, we support [GitHub Codespaces](https://github.com/codespaces): [Create a new Codespace](https://github.com/codespaces/new?machine=basicLinux32gb&repo=334825694&ref=develop&devcontainer_path=.devcontainer%2Fdevcontainer.json&location=SouthEastAsia)

We also support [Gitpod](https://www.gitpod.io/): [Open in Gitpod](https://gitpod.io/#https://github.com/deepmodeling/abacus-develop)

## Install by conda

Conda is a package management system with a separated environment, not requiring system privileges. A pre-built ABACUS binary with all requirements is available at [conda-forge](https://anaconda.org/conda-forge/abacus). It supports advanced features including Libxc, LibRI, and DeePKS. Conda will install the GPU-supported version of ABACUS if a valid GPU driver is present. Please refer to [the advanced installation guide](../advanced/install.md) for more details.

```bash
# Install
# We recommend installing ABACUS in a new environment to avoid potential conflicts:
conda create -n abacus_env abacus -c conda-forge

# Run
conda activate abacus_env
OMP_NUM_THREADS=1 mpirun -n 4 abacus

# Update
conda update -n abacus_env abacus -c conda-forge
```

> If OpenBLAS gives warning about OpenMP threads, please install conda package `openblas=*=openmp*` or `blas=*=mkl`. See [switching BLAS implementation in conda](https://conda-forge.org/docs/maintainer/knowledge_base.html#switching-blas-implementation).

> ABACUS supports `OpenMPI` and `MPICH` variant. Install `mpich` or `openmpi` package to switch MPI library if required.

For more details on building a conda package of ABACUS locally, please refer to the [conda recipe file](https://github.com/deepmodeling/abacus-develop/blob/develop/conda/meta.yaml).

> Note: The [deepmodeling conda channel](https://anaconda.org/deepmodeling/abacus) offers historical versions of ABACUS.

### Developing with conda

It is possible to build ABACUS from source based on the conda environment.

```bash
conda create -n abacus_env abacus -c conda-forge
conda activate abacus_env
export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH

# By default OpenBLAS is used; run `conda install "blas=*=mkl" mkl_fft mkl-devel -c conda-forge` to switch implementation.
export MKLROOT=$CONDA_PREFIX # If Intel MKL is required.

export CMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'`:$CMAKE_PREFIX_PATH # If DEEPKS support is required;
# usually expands to `$CONDA_PREFIX/lib/python3.1/site-packages/torch/share/cmake`
```

And, follow the instructions in [Build and Install](#build-and-install) part above withou manually setting paths to dependencies.
See [the advanced installation guide](../advanced/install.md) for more features.
Make sure the environment variables are set before running `cmake`.
Possible command: `cmake -B build -DENABLE_DEEPKS=ON -DENABLE_LIBXC=ON -DENABLE_LIBRI=ON`.
