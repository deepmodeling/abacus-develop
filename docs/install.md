# Download and install

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Build and install ABACUS with CMake](#build-and-install-abacus-with-cmake)
  - [Build ABACUS with make](#build-abacus-with-make)
    - [Link LIBXC](#link-libxc)
- [Installation with DeePKS](#installation-with-deepks)
  - [Extra prerequisites](#extra-prerequisites)
  - [Extra settings for building](#extra-settings-for-building)

  [back to main page](../README.md)

## Installation

### Container Deployment

We've built a ready-for-use version of ABACUS with docker [here](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus). For a quick start: pull the image, prepare the data, run container. Instructions on using the image can be accessed in [Dockerfile](../Dockerfile). A mirror is available by `docker pull registry.dp.tech/deepmodeling/abacus`.

We also offer a pre-built docker image containing all the requirements for development. Please refer to our [Package Page](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus-development-kit).

The project is ready for VS Code development container. Please refer to [Developing inside a Container](https://code.visualstudio.com/docs/remote/containers#_quick-start-try-a-development-container). Choose `Open a Remote Window -> Clone a Repository in Container Volume` in VS Code command palette, and put the [git address](https://github.com/deepmodeling/abacus-develop.git) of `ABACUS` when prompted.

We also support [Gitpod](https://www.gitpod.io/) to offer an ready-to-use online development environment.
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/deepmodeling/abacus-develop)

Please note that containers target at developing and testing, but not massively parallel computing. Docker has a bad support to MPI; for production, please compile ABACUS from source code to avoid compatibility issues and gain a better performace.

### Prerequisites

ABACUS currently supports Linux. `Dockerfile`s under the root directory of our repo will come in handy.

To compile ABACUS, please make sure that the following prerequisites are present:

- C++ compiler, supporting C++11. You can use [Intel® C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/).
- MPI compiler. The recommended version are [Intel MPI](https://software.intel.com/enus/mpi-library) or [MPICH](https://www.mpich.org/).
- Fortran compiler if you are building `BLAS`, `LAPACK`, `ScaLAPACK`, and `ELPA` from source file. You can use [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or [GFortran](https://gcc.gnu.org/fortran/).
- [BLAS](http://www.netlib.org/blas/). You can use [OpenBLAS](https://www.openblas.net/).
- [LAPACK](http://www.netlib.org/lapack/).
- [ScaLAPACK](http://www.netlib.org/scalapack/).
- [FFTW3](http://www.fftw.org/).
- [ELPA](https://elpa.mpcdf.mpg.de/) >= 2017.
- [CEREAL](https://uscilab.github.io/cereal/).
- [Libxc](https://tddft.org/programs/libxc/)>=5.1.7(optional), if you wish to use meta-GGA functionals.

> GCC version 5 or later is required; Intel compilers also use GCC headers and libraries[(ref)](https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compatibility-and-portability/gcc-compatibility-and-interoperability.html#gcc-compatibility-and-interoperability_GUID-52CB6FE0-83DA-4028-9EF4-0DFAF1652736).

These packages can be installed with popular package management system, such as `apt` and `yum`:

```bash
sudo apt update && sudo apt install -y libopenblas-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc git
```

> Installing ELPA by apt only matches requirements on Ubuntu 22.04. For earlier linux distributions, you should install elpa from source.

Alternatively, you can choose [Intel® oneAPI toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/commercial-base-hpc.html) (former Parallel Studio) as toolchain. The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#base-kit) contains Intel® oneAPI Math Kernel Library (aka `MKL`), including `BLAS`, `LAPACK`, `ScaLAPACK` and `FFTW3`. The [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#hpc-kit) contains Intel® MPI Library, and C++ compiler(including MPI compiler). Please noted that building `elpa` with a different MPI library may cause conflict between MPI libraries. Don't forget to [set environment variables](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-render-linux/top/configure-your-system.html) before you start! `cmake` will use Intel MKL if the environment variable `MKLROOT` is set.

> Please refer to our [guide](https://github.com/deepmodeling/abacus-develop/wiki/Building-and-Running-ABACUS) on requirements.

If you have trouble building requirements, our Dockerfiles in root path offer a reference, or read the section below to use a pre-built container.

And of course, a copy of ABACUS source code is required:

- Clone the whole repo with git: `git clone https://github.com/deepmodeling/abacus-develop.git`
- Clone the minimum required part of repo: `git clone https://github.com/deepmodeling/abacus-develop.git --depth=1`
- Download the latest source code without git: `wget https://github.com/deepmodeling/abacus-develop/archive/refs/heads/develop.zip`
- Get the source code of a stable version [here](https://github.com/deepmodeling/abacus-develop/releases)
- If you have connection issues accessing GitHub, please try out our official [Gitee repo](https://gitee.com/deepmodeling/abacus-develop/): replacing 'github.com' with 'gitee.com' works for all the links above. e.g. `git clone https://gitee.com/deepmodeling/abacus-develop.git`


[back to top](#download-and-install)

### Build and install ABACUS with CMake

We recommend building ABACUS with `cmake` to avoid dependency issues. `Makefile` is deprecated.

#### Configure

ABACUS requires a minimum `cmake` version of `3.16`, and `3.18` with advanced features like the integration with DeePKS or utilizing GPU. Check the version of `cmake` on your machine with:

```bash
cmake --version
```

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`. If no install prefix is specified, the binary will be installed to `/usr/local/bin/abacus` by default.

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH}
```

You can provide path of each dependent package if the package cannot be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDE_DIR`, `MPI_CXX_COMPILER` and `MKLROOT` are currently available to specify.
For example:

```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3
```

If environment variable `MKLROOT` exists, `cmake` will take MKL as a preference, i.e. not using `LAPACK` and `ScaLAPACK`. To disable MKL, unset environment variable `MKLROOT`, or pass `-DMKLROOT=OFF` to `cmake`.

You can also choose to build with which components, e.g.:

```bash
cmake -B build -DENABLE_LIBXC=1 -DUSE_CUDA=1
```

If Libxc is not installed in standard path (i.e. installed with a custom prefix path), you can set `Libxc_DIR` to the corresponding directory.

```bash
cmake -B build -DLibxc_DIR=~/libxc
```

To build tests for ABACUS, define `BUILD_TESTING` flag. You can also specify path to local installation of [Googletest](https://github.com/google/googletest) by setting `GTEST_DIR` flags. If not found in local, the configuration process will try to download it automatically.

```bash
cmake -B build -DBUILD_TESTING=1
```

#### Build and Install

After configuring, start build and install by:

```bash
cmake --build build -j`nproc`
cmake --install build
```

You can change the number after `-j` on your need: set to the number of CPU cores(`nproc`) to gain the best performance.

[back to top](#download-and-install)

### Build ABACUS with make

<!-- Before starting to build the program, note that if you are using Intel MKL library, please set the following environmental variable:

```bash
export MKL_NUM_THREAD=1
``` -->

> Note: Makefile is not always updated in time. We suggest using CMake to configure and compile.

To compile the ABACUS program using legacy `make`, uses only need to edit the file `Makefile.vars` under `source` directory:

```bash
cd source/
vi Makefile.vars
```

Specify the location of the compiler and libraries present in your own machine:

```makefile
# This is the Makefile of ABACUS API
#======================================================================
# Users set
#======================================================================
CC = mpiicpc
# mpiicpc:   compile intel parallel version
# icpc:      compile gnu serial version
# make: ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: nothing need to be set except LIBXC_DIR
# 
# mpicxx:    compile gnu parallel version
# g++:       compile gnu serial version
# make: FFTW_DIR, OPENBLAS_LIB_DIR, SCALAPACK_LIB_DIR, ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: FFTW_DIR, OPENBLAS_LIB_DIR must be set.
#======================================================================

#-------  FOR INTEL COMPILER  ------------
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
# directory of elpa, which contains include and lib/libelpa.a

CEREAL_DIR    = /public/soft/cereal
# directory of cereal, which contains a include directory in it.

#-------  FOR GNU COMPILER  ---------------
# FFTW_DIR = /public/soft/fftw_3.3.8
# # directory of fftw package, which contains lib/libfftw3.a. Only used when CC = mpicxx/g++

# OPENBLAS_LIB_DIR   = /public/soft/openblas/lib
# # directory of libopenblas.a, only used when CC = mpicxx/g++

# SCALAPACK_LIB_DIR  = /public/soft/openblas/lib
# # directory of libscalapack.a, only used when CC = mpicxx/g++

# ELPA_DIR      = /public/soft/elpa_21.05.002
# ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
# # directory of elpa, which contains include and lib/libelpa.a

# CEREAL_DIR    = /public/soft/cereal
# # directory of cereal, which contains a include directory in it.

#------  OPTIONAL LIBS  -----------

# LIBTORCH_DIR  = /usr/local
# LIBNPY_DIR    = /usr/local
# add them to use DEEPKS

# LIBXC_DIR    		= /public/soft/libxc
# directory of libxc(>5.1.7), which contains include and lib/libxc.a
# add LIBXC_DIR to use libxc to compile ABACUS
#======================================================================
```

For example, below is a case where the Intel C++ compiler, Intel MPI and CEREAL_DIR are used, along with Intel MKL library. The file Makefile.vars can be set as
follows:

```makefile
CC = mpiicpc/icpc
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR    = /public/soft/cereal
```
When `CC=mpiicpc`, a parallel version will be compiled. When `CC=icpc`, a serial version will be compiled.


Another example is where GCC, GFORTRAN, MPICH, ScaLAPACK, ELPA and CEREAL are used:

```makefile
CC = mpicxx/g++
FFTW_DIR = /public/soft/fftw_3.3.8
OPENBLAS_LIB_DIR   = /public/soft/openblas/lib
SCALAPACK_LIB_DIR  = /public/soft/openblas/lib
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR    = /public/soft/cereal
```
When `CC=mpicxx`, a parallel version will be compiled. When `CC=g++`, a serial version will be compiled.

Except modifying `Makefile.vars`, you can also directly use
```makefile
make CC=mpicpc ELPA_DIR=/public/soft/elpa_21.05.002 \
ELPA_INCLUDE_DIR=${ELPA_DIR}/include/elpa-2021.05.002 \
CEREAL_DIR=/public/soft/cereal
```
ABACUS now support full version and pw version. Use `make` or `make abacus` to compile full version which supports LCAO calculations. Use `make pw` to compile pw version which only supports pw calculations. For pw version, `make pw CC=mpiicpc`, you do not need to provide any libs. For `make pw CC=mpicxx`, you need provide `FFTW_DIR` and `OPENBLAS_LIB_DIR`.

Besides, libxc and deepks are optional libs to compile abacus. 
They will be used when `LIBXC_DIR` is defined like
```
LIBXC_DIR    		= /public/soft/libxc
```
or `LIBTORCH_DIR` and `LIBNPY_DIR` like
```makefile
LIBTORCH_DIR  = /usr/local
LIBNPY_DIR    = /usr/local
```

After modifying the `Makefile.vars` file, execute `make` or `make -j12` or `make -j`to build the program.

After the compilation finishes without error messages (except perhaps for some warnings), an executable program `ABACUS.mpi` will be created in directory `bin/`.

[back to top](#download-and-install)

#### Link LIBXC

The program compiled using the above instructions do not link with LIBXC and use exchange-correlation functionals as written in the ABACUS program. However, for some functionals (such as HSE hybrid functional), LIBXC is required.

To compile ABACUS with LIBXC, you need to define `LIBXC_DIR` in the file `Makefile.vars` or use 
```makefile
make LIBXC_DIR=/pulic/soft/libxc
``` 
directly.

[back to top](#download-and-install)

## Installation with DeePKS

This part of installation is based on [Installation](#installation). If DeePKS feature is requied for [DeePKS-kit](https://github.com/deepmodeling/deepks-kit), the following prerequisites and steps are needed:

### Extra prerequisites

- C++ compiler, supporting **C++14**
- [LibTorch](https://pytorch.org/) with cxx11 ABI supporting CPU
- [Libnpy](https://github.com/llohse/libnpy/)

### Extra settings for building

### Using Cmake

```bash
cmake -B build -DENABLE_DEEPKS=1
```

### Using Makefile

To compile ABACUS with DEEPKS, you need to define `LIBTORCH_DIR` and `LIBNPY_DIR` in the file `Makefile.vars` or use 
```makefile
make LIBTORCH_DIR=/opt/libtorch/ LIBNPY_DIR=/opt/libnpy/
``` 
directly.


[back to top](#download-and-install)
