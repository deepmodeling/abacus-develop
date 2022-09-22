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

### Add Libxc Support
To compile ABACUS with Libxc, modifications should be made in three files:

First of all, in the file `Makefile.vars`, apart from the variables above, further provide the location of Libxc:

```bash
LIBXC_DIR =
```

Then, in the file 'Makefile.system', add "${Libxc_LIB}" to the `LIBS` flag, for example:

```bash
LIBS = -lifcore -lm -lpthread ${LAPACK_LIB} ${FFTW_LIB} ${ELPA_LIB} ${Libxc_LIB}
```

Finally, in `Makefile`, add "-DUSE_Libxc" to the `HONG` flag, for example:

```bash
HONG_MPI_SELINV_20210523 = -D__FP ${HONG_FFTW} ${HONG_LAPACK} -D__LCAO -D__MPI -D__OPENMP -D__SELINV -DMETIS -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1 -DUSE_CEREAL_SERIALIZATION -D__EXX -DUSE_Libxc
HONG=${HONG_MPI_SELINV_20210523}
```
### Add DeePKS Support


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

In `Makefile`, set the Macro as `HONG_DEEPKS`:

```Makefile
#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!
#! change series version or parallel version~~~
#HONG=${HONG_MPI_SELINV_20210523}
#HONG=${HONG_SER_SELINV}
HONG=${HONG_DEEPKS}
```
