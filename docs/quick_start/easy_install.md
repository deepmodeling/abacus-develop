# Easy Installation

This guide helps you install ABACUS with basic features. For DeePKS, Libxc support or building with `make`, please refer to [the advanced installation guide](../advanced/install.md). We recommend building ABACUS with `cmake` to avoid dependency issues.

## Prerequisites

ABACUS currently supports Linux. To compile ABACUS, please make sure that the following prerequisites are present:

- C++ compiler, supporting C++11. You can use [Intel® C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/).
- MPI compiler. The recommended version are [Intel MPI](https://software.intel.com/enus/mpi-library) or [MPICH](https://www.mpich.org/).
- Fortran compiler if you are building `BLAS`, `LAPACK`, `ScaLAPACK`, and `ELPA` from source file. You can use [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or [GFortran](https://gcc.gnu.org/fortran/).
- [BLAS](http://www.netlib.org/blas/). You can use [OpenBLAS](https://www.openblas.net/).
- [LAPACK](http://www.netlib.org/lapack/).
- [ScaLAPACK](http://www.netlib.org/scalapack/).
- [FFTW3](http://www.fftw.org/).
- [ELPA](https://elpa.mpcdf.mpg.de/) >= 2017.
- [CEREAL](https://uscilab.github.io/cereal/).

> GCC version 5 or later is required; Intel compilers also use GCC headers and libraries[(ref)](https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compatibility-and-portability/gcc-compatibility-and-interoperability.html#gcc-compatibility-and-interoperability_GUID-52CB6FE0-83DA-4028-9EF4-0DFAF1652736).

These packages can be installed with popular package management system, such as `apt` and `yum`:

```bash
sudo apt update && sudo apt install -y libopenblas-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc git
```

> Installing ELPA by apt only matches requirements on Ubuntu 22.04. For earlier linux distributions, you should install elpa from source.

Alternatively, you can choose [Intel® oneAPI toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/commercial-base-hpc.html) (former Parallel Studio) as toolchain. The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#base-kit) contains Intel® oneAPI Math Kernel Library (aka `MKL`), including `BLAS`, `LAPACK`, `ScaLAPACK` and `FFTW3`. The [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#hpc-kit) contains Intel® MPI Library, and C++ compiler(including MPI compiler). Please noted that building `elpa` with a different MPI library may cause conflict between MPI libraries. Don't forget to [set environment variables](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-render-linux/top/configure-your-system.html) before you start! `cmake` will use Intel MKL if the environment variable `MKLROOT` is set.

> Please refer to our [guide](https://github.com/deepmodeling/abacus-develop/wiki/Building-and-Running-ABACUS) on requirements.

And of course, a copy of ABACUS source code is required:

- Clone the whole repo with git: `git clone https://github.com/deepmodeling/abacus-develop.git`
- Clone the minimum required part of repo: `git clone https://github.com/deepmodeling/abacus-develop.git --depth=1`
- Download the latest source code without git: `wget https://github.com/deepmodeling/abacus-develop/archive/refs/heads/develop.zip`
- Get the source code of a stable version [here](https://github.com/deepmodeling/abacus-develop/releases)
- If you have connection issues accessing GitHub, please try out our official [Gitee repo](https://gitee.com/deepmodeling/abacus-develop/): replacing 'github.com' with 'gitee.com' works for all the links above. e.g. `git clone https://gitee.com/deepmodeling/abacus-develop.git`

## Configure

ABACUS requires a minimum `cmake` version of `3.16`. Check the version of `cmake` on your machine with:

```bash
cmake --version
```

You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`. If no install prefix is specified, the binary will be installed to `/usr/local/bin/abacus` by default.

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=${YOUR_PATH_TO_ABACUS_BINARY}
```

You can provide path of each dependent package if the package can not be automatically found by cmake.
Keys `LAPACK_DIR`, `SCALAPACK_DIR`, `ELPA_DIR`, `FFTW3_DIR`, `CEREAL_INCLUDE_DIR`, `MPI_CXX_COMPILER` and `MKLROOT` are currently available to specify.
For example:

```bash
cmake -B build -DFFTW3_ROOT=/opt/fftw3
```

If environment variable `MKLROOT` exists, `cmake` will take MKL as a preference, i.e. not using `LAPACK` and `ScaLAPACK`. To disable MKL, unset environment variable `MKLROOT`, or pass `-DMKLROOT=OFF` to `cmake`.


## Build and Install

After configuring, start build and install by:

```bash
cmake --build build -j`nproc`
cmake --install build
```

You can change the number after `-j` on your need: set to the number of CPU cores(`nproc`) to reduce compilation time.

## Container Deployment

We've built a ready-for-use version of ABACUS with docker [here](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus). For a quick start: pull the image, prepare the data, run container. Instructions on using the image can be accessed in [Dockerfile](../../Dockerfile). A mirror is available by `docker pull registry.dp.tech/deepmodeling/abacus`.

We also offer a pre-built docker image containing all the requirements for development. Please refer to our [Package Page](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus-development-kit).

The project is ready for VS Code development container. Please refer to [Developing inside a Container](https://code.visualstudio.com/docs/remote/containers#_quick-start-try-a-development-container). Choose `Open a Remote Window -> Clone a Repository in Container Volume` in VS Code command palette, and put the [git address](https://github.com/deepmodeling/abacus-develop.git) of `ABACUS` when prompted.

We also support [Gitpod](https://www.gitpod.io/) to offer an ready-to-use online development environment.
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/deepmodeling/abacus-develop)

> Please note that containers target at developing and testing, but not massively parallel computing. Docker has a bad support to MPI; for production, please compile ABACUS from source code to avoid compatibility issues and gain a better performace.
