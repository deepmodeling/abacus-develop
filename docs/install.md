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

### Prerequisites

ABACUS current supports Linux. To compile ABACUS, please make sure that the following prerequisites are
present:

- C++ compiler, supporting C++11. You can use [Intel® C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/).
- MPI compiler. The recommended version are [Intel MPI](https://software.intel.com/enus/mpi-library) or [MPICH](https://www.mpich.org/).
- Fortran compiler for building `BLAS`, `LAPACK` and `ScaLAPACK`. You can use[Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) [GFortran](https://gcc.gnu.org/fortran/).
- [BLAS](http://www.netlib.org/blas/). You can use [OpenBLAS](https://www.openblas.net/).
- [LAPACK](http://www.netlib.org/lapack/).
- [ScaLAPACK](http://www.netlib.org/scalapack/).
- [FFTW](http://www.fftw.org/).
- [ELPA](https://elpa.mpcdf.mpg.de/).
- [CEREAL](https://uscilab.github.io/cereal/).

Alternatively, you can choose [Intel® oneAPI toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/commercial-base-hpc.html) (former Parallel Studio) as toolchain. The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#base-kit) contains Intel® oneAPI Math Kernel Library (aka `MKL`), including `BLAS`, `LAPACK`, `ScaLAPACK` and `FFTW3`,  - this means that no Fortran compiler required anymore. The [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#hpc-kit) contains Intel® MPI Library, and C++ compiler(including MPI compiler). Please noted that building `elpa` with a different MPI library may cause conflict between MPI libraries. Don't forget to [set environment variables](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-render-linux/top/configure-your-system.html) before you start! `cmake` will use Intel MKL if the environment variable `MKLROOT` is set.

If you have trouble building requirements, our Dockerfiles in root path offer a reference, or read the section below to use a pre-built container.
[back to top](#download-and-install)

### Container Deployment

We offer a pre-built docker image containing all the requirements - you only need to clone and compile `abacus` in the container. Please refer to our [Package Page](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus-development-kit).

The project is ready for VS Code development container. Please refer to [Developing inside a Container](https://code.visualstudio.com/docs/remote/containers#_quick-start-try-a-development-container). Choose `Open a Remote Window -> Clone a Repository in Container Volume` in VS Code command palette, and put the [git address](https://github.com/deepmodeling/abacus-develop.git) of `ABACUS` when prompted.

We also support [gitpod](https://www.gitpod.io/) to offer an ready-to-use online development environment.
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/deepmodeling/abacus-develop)

### Build and install ABACUS with CMake

We recommend building ABACUS with `cmake` to avoid dependency issues. `Makefile` is deprecated.

#### Configure

ABACUS requires a minimum `cmake` version of `3.18`. Check the version of `cmake`  on your machine with:

```bash
cmake --version
```

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`. If no install prefix is specified, the binary will be installed to `/usr/local/bin/abacus` by default.

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH}
```

You can provide path of each dependent package if the package cannot be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDEDIR`, `MPI_CXX_COMPILER` and `MKLROOT` are currently available to specify.
For example:

```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3
```

If environment variable `MKLROOT` exists, `cmake` will take MKL as a preference, i.e. not using `LAPACK` and `ScaLAPACK`. To disable MKL, unset environment variable `MKLROOT`, or pass `-DMKLROOT=OFF` to `cmake`.

You can also choose to build with which components.

```bash
cmake -B build -DUSE_LIBXC=1 -DUSE_CUDA=1
```

If Libxc is not installed in standard path (i.e. installed with a custom prefix path), you may add the installation prefix of `FindLibxc.cmake` to `CMAKE_MODULE_PATH` environment variable, or set `Libxc_DIR` to the directory containing the file.

```bash
cmake -B build -DLibxc_DIR=~/libxc
```

To build tests for abacus, define `BUILD_TESTING` flag. You can also specify path to local installation of [Googletest](https://github.com/google/googletest) by setting `GTEST_DIR` flags. If not found in local, the configuration process will try to download it automatically.

```bash
cmake -B build -DBUILD_TESTING=1
```

#### Build and Install

After configuring, start build and install by:

```bash
cmake --build build -j9
cmake --install build
```

`-j9` specifies the number of jobs to run simultaneously. You can change the number on your need: set to the number of CPU cores to gain the best performance.

[back to top](#download-and-install)

### Build ABACUS with make

<!-- Before starting to build the program, note that if you are using Intel MKL library, please set the following environmental variable:

```bash
export MKL_NUM_THREAD=1
``` -->

To compile the ABACUS program using legacy `make`, first edit the file `Makefile.vars` under `source` directory:

```bash
cd source/
vi Makefile.vars
```

Specify the location of the compiler and libraries present in your own machine:

```bash
CPLUSPLUS =
CPLUSPLUS_MPI =
FORTRAN =
LAPACK_DIR =
FFTW_DIR =
BOOST_DIR =
ELPA_DIR =
CEREAL_DIR =
```

For example, below is a case where the Intel C++ compiler, Intel MPI are used, along with Intel MKL library. The file Makefile.vars can be set as
follows:

```bash
CPLUSPLUS = icpc
CPLUSPLUS_MPI = mpiicpc
FORTRAN = ifort
LAPACK_DIR = /opt/intel/.../mkl/lib/intel64/
FFTW_DIR = /opt/fftw-3.3.8/
BOOST_DIR = /opt/boost/1.64.0/
ELPA_DIR = /opt/elpa/2016.05.004/
CEREAL_DIR = /opt/cereal/
```

Another example is where GCC, GFORTRAN, MPICH and ScaLAPACK are used:

```bash
CPLUSPLUS = g++
CPLUSPLUS_MPI = mpicxx
FORTRAN = gfortran
SCALAPACK_DIR = /opt/scalapack/
FFTW3_DIR = /opt/fftw-3.3.8/
BOOST_DIR = /opt/boost/1.64.0/
ELPA_DIR = /opt/elpa/2016.05.004/
CEREAL_DIR = /opt/cereal/
```

For this option, it is further required to set the parameter `LIBS` in `Makefile.system`:

```bash
LIBS = \
  -lgfortran -lm \
  -openmp -lpthread \
  ${SCALAPACK_DIR}/lib/libscalapack.a \
  /opt/lapack/lib/liblapack.a \
  /opt/blas/lib/libblas.a \
  /opt/blacs/lib/libblacs.a \
  ${FFTW_LIB} \
  ${ELPA_LIB} \

```

After modifying the `Makefile.vars` file, execute `make` to build the program.

```bash
make -j
```

After the compilation finishes without error messages (except perhaps for some warnings), an executable program `ABACUS.mpi` will be created in directory `bin/`.

[back to top](#download-and-install)

#### Link LIBXC

The program compiled using the above instructions do not link with LIBXC and use exchange-correlation functionals as written in the ABACUS program. However, for some functionals (such as HSE hybrid functional), LIBXC is required.

To compile ABACUS with LIBXC, modifications should be made in three files:

First of all, in the file `Makefile.vars`, apart from the variables above, further provide the location of LIBXC:

```bash
LIBXC_DIR =
```

Then, in the file 'Makefile.system', add "${LIBXC_LIB}" to the `LIBS` flag, for example:

```bash
LIBS = -lifcore -lm -lpthread ${LAPACK_LIB} ${FFTW_LIB} ${ELPA_LIB} ${LIBXC_LIB}
```

Finally, in `Makefile`, add "-DUSE_LIBXC" to the `HONG` flag, for example:

```bash
HONG_MPI_SELINV_20210523 = -D__FP ${HONG_FFTW} ${HONG_LAPACK} -D__LCAO -D__MPI -D__OPENMP -D__SELINV -DMETIS -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1 -DUSE_CEREAL_SERIALIZATION -D__EXX -DUSE_LIBXC
HONG=${HONG_MPI_SELINV_20210523}
```

[back to top](#download-and-install)

## Installation with DeePKS

This part of installation is based on [Installation](#installation). If DeePKS feature is requied for [DeePKS-kit](https://github.com/deepmodeling/deepks-kit), the following prerequisites and steps are needed:

### Extra prerequisites

- C++ compiler, supporting **C++14**. For example, Intel C++ compiler 18
- [LibTorch](https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcpu.zip) for cpu, with c++11 ABI;
- [Libnpy](https://github.com.cnpmjs.org/llohse/libnpy/);

### Extra settings for building

### Using Cmake

```bash
cmake -B build -DENABLE_DEEPKS=1
```

### Using Makefile

Set `LIBTORCH_DIR`and `LIBNPY_DIR`in `Makefile.vars`. For example:

```Makefile
LIBTORCH_DIR = /opt/libtorch/
LIBNPY_DIR = /opt/libnpy/
```

In `Makefile.system`, add `LIBTORCH_LIB` to  `LIBS`, then set `-std=c++14` in `OPTS`:

```Makefile
LIBS = -lifcore -lm -lpthread ${LIBTORCH_LIB} ${LAPACK_LIB} ${FFTW_LIB} ${ELPA_LIB} # for DeePKS
#LIBS = -lifcore -lm -lpthread ${LAPACK_LIB} ${FFTW_LIB} ${ELPA_LIB}
```

```Makefile
OPTS = ${INCLUDES} -Ofast -traceback -std=c++14 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g
```

- module_base
- module_cell
- module_grid
- module_md
- module_neighbor
- module_orbital
- obj
- src_external
- src_global
- src_io
- src_ions
- src_lcao
- src_parallel
- src_pdiag
- src_pw
- src_ri

In `Makefile`, set the Macro as `HONG_DEEPKS`:

```Makefile
#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!
#! change series version or parallel version~~~
#HONG=${HONG_MPI_SELINV_20210523}
#HONG=${HONG_SER_SELINV}
HONG=${HONG_DEEPKS}
```

[back to top](#download-and-install)
