# Advanced Installation Options
This guide helps you install ABACUS with advanced features. Please make sure to read the [easy-installation guide](../quick_start/easy_install.md) before.
## Build with Libxc

ABACUS use exchange-correlation functionals by default. However, for some functionals (such as HSE hybrid functional), Libxc is required.

Dependency: [Libxc](https://tddft.org/programs/libxc/)>=5.1.7 .

If Libxc is not installed in standard path (i.e. installed with a custom prefix path), you can set `LIBXC_DIR` to the corresponding directory.

```bash
cmake -B build -DLIBXC_DIR=~/libxc
```

## Build with DeePKS
If DeePKS feature is requied for [DeePKS-kit](https://github.com/deepmodeling/deepks-kit), the following prerequisites and steps are needed:

### Extra prerequisites

- C++ compiler, supporting **C++14**
- CMake `3.18` and above
- [LibTorch](https://pytorch.org/) with cxx11 ABI supporting CPU
- [Libnpy](https://github.com/llohse/libnpy/)

```bash
cmake -B build -DENABLE_DEEPKS=1 -DTorch_DIR=~/libtorch/share/cmake/Torch/ -Dlibnpy_INCLUDE_DIR=~/libnpy/include
```

## Build Unit Tests
To build tests for ABACUS, define `BUILD_TESTING` flag. You can also specify path to local installation of [Googletest](https://github.com/google/googletest) by setting `GTEST_DIR` flags. If not found in local, the configuration process will try to download it automatically.

```bash
cmake -B build -DBUILD_TESTING=1
```
## Build ABACUS with make

> Note: We suggest using CMake to configure and compile.

To compile the ABACUS program using legacy `make`, users only need to edit the file `Makefile.vars` under `source` directory:

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
# icpc:      compile intel serial version
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

For example, below is a case where the Intel C++ compiler, Intel MPI and CEREAL are used, along with Intel MKL library. The file Makefile.vars can be set as
follows:

```makefile
CC = mpiicpc #(or CC = icpc)
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR    = /public/soft/cereal
```
When `CC=mpiicpc`, a parallel version will be compiled. When `CC=icpc`, a serial version will be compiled.


Another example is where the Gnu C++ compiler, MPICH, OPENBLAS, ScaLAPACK, ELPA and CEREAL are used:

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
make CC=mpiicpc ELPA_DIR=/public/soft/elpa_21.05.002 \
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

### Add Libxc Support

The program compiled using the above instructions do not link with LIBXC and use exchange-correlation functionals as written in the ABACUS program. However, for some functionals (such as HSE hybrid functional), LIBXC is required.

To compile ABACUS with LIBXC, you need to define `LIBXC_DIR` in the file `Makefile.vars` or use 
```makefile
make LIBXC_DIR=/pulic/soft/libxc
``` 
directly.

### Add DeePKS Support

To compile ABACUS with DEEPKS, you need to define `LIBTORCH_DIR` and `LIBNPY_DIR` in the file `Makefile.vars` or use 
```makefile
make LIBTORCH_DIR=/opt/libtorch/ LIBNPY_DIR=/opt/libnpy/
``` 
directly.
