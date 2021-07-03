# Download and install

- [Structure of the package](#structure-of-the-package)
  - [Structure of source code](#structure-of-source-code)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Building the program](#building-the-program)

  [back to main page](#../install.md)  

# Structure of the package
Under the ABACUS directory, there are the following subdirectories:

- cmake/

  which contains relevant files for compiling the code with cmake
- documents/

  which contains a copy of the manual in pdf format
- examples/

  which contains some examples
- source/

  which contains the source code and makefiles
- tests/

  which contains test examples
- tools/

  which currently contains the script for generating the numerical atomic orbitals

[back to top](#download-and-install)

## Structure of source code
The source directory further contains the following folders, where the source files of ABACUS are located:
- module_base
- module_cell
- module_grid
- module_grid
- module_neighbor
- module_ORB
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

[back to top](#download-and-install)

# Installation

## Prerequisites
In order to compile ABACUS, users should make sure that the following prerequisites are
present:

- C++ compiler, supporting C++11. For example, [Intel C++ compiler](https://software.intel.com/enus/c-compilers) or [GCC](https://gcc.gnu.org/);
- Fortran compiler;
- MPI compiler. The recommended version are [Intel MPI](https://software.intel.com/enus/mpi-library) or [MPICH](https://www.mpich.org/);
- [Boost C++ library](https://www.boost.org/);
- The ScaLAPACK library. For example, [Intel MKL](https://software.intel.com/en-us/mkl)
or [Netlib ScaLAPACK](http://www.netlib.org/scalapack/);
- The [FFTW library](http://www.fftw.org/). ABACUS now supports both FFTW2 and
FFTW3;
- The [ELPA library](https://elpa.mpcdf.mpg.de/);
- The [CEREAL library](https://uscilab.github.io/cereal/);

[back to top](#download-and-install)

## Building the program
To compile the ABACUS program, go to the source directory:
```bash
cd source/
```
Then open and edit the file `Makefile.vars` using any editor tools you like, e.g., vi:
```bash
vi Makefile.vars
```
Specify the location of the compiler and libraries present in your own machine:
```
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
```
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
```
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
```
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
After modifying the `Makefile.vars` file, to build the program, simply type
```bash
make
```
If the compilation finishes without error messages (except perhaps for some warnings), an executable program `ABACUS.mpi` will be created in directory `bin/`

[back to top](#download-and-install)