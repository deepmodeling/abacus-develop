# Advanced Installation Options

This guide helps you install ABACUS with advanced features. Please make sure to read the [easy-installation guide](../quick_start/easy_install.md) before.

## Build with Libxc

ABACUS use exchange-correlation functionals by default. However, for some functionals (such as HSE hybrid functional), Libxc is required.

Dependency: [Libxc](https://tddft.org/programs/libxc/) >= 5.1.7 .

> Note: Building Libxc from source with Makefile does NOT support using it in CMake here. Please compile Libxc with CMake instead.

If Libxc is not installed in standard path (i.e. installed with a custom prefix path), you can set `Libxc_DIR` to the corresponding directory.

```bash
cmake -B build -DLibxc_DIR=~/libxc
```

## Build with DeePKS

If DeePKS feature is required for [DeePKS-kit](https://github.com/deepmodeling/deepks-kit), the following prerequisites and steps are needed:

- C++ compiler, supporting **C++14** (GCC >= 5 is sufficient)
- CMake >= 3.18
- [LibTorch](https://pytorch.org/) with cxx11 ABI supporting CPU
- [Libnpy](https://github.com/llohse/libnpy/)

```bash
cmake -B build -DENABLE_DEEPKS=1 -DTorch_DIR=~/libtorch/share/cmake/Torch/ -Dlibnpy_INCLUDE_DIR=~/libnpy/include
```

> CMake will try to download Libnpy if it cannot be found locally.

## Build with DeePMD-kit

> Note: This part is only required if you want to load a trained DeeP Potential and run molecular dynamics with that. To train the DeeP Potential with DP-GEN, no extra prerequisite is needed and please refer to [this page](http://abacus.deepmodeling.com/en/latest/advanced/interface/dpgen.html) for ABACUS interface with DP-GEN.

If the Deep Potential model is employed in Molecule Dynamics calculations, the following prerequisites and steps are needed:

- [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit)
- [TensorFlow](https://www.tensorflow.org/)

```bash
cmake -B build -DDeePMD_DIR=~/deepmd-kit -DTensorFlow_DIR=~/tensorflow
```

> `deepmd_c`/`deepmd_cc` and `tensorflow_cc` libraries would be called according to `DeePMD_DIR` and `TensorFlow_DIR`, which is showed in detail in [this page](https://github.com/deepmodeling/deepmd-kit/blob/master/doc/inference/cxx.md).

## Build with LibRI and LibComm

The new EXX implementation depends on two external libraries:

- [LibRI](https://github.com/abacusmodeling/LibRI)
- [LibComm](https://github.com/abacusmodeling/LibComm)

These two libraries are added as submodules in the [deps](https://github.com/deepmodeling/abacus-develop/tree/develop/deps) folder. Set `-DENABLE_LIBRI=ON` to build with these two libraries.

If you prefer using manually downloaded libraries, provide `-DLIBRI_DIR=${path to your LibRI folder} -DLIBCOMM_DIR=${path to your LibComm folder}`. 

## Build Unit Tests

To build tests for ABACUS, define `BUILD_TESTING` flag. You can also specify path to local installation of [Googletest](https://github.com/google/googletest) by setting `GTEST_DIR` flags. If not found in local, the configuration process will try to download it automatically.

```bash
cmake -B build -DBUILD_TESTING=1
```

After building and installing, unit tests can be performed with `ctest`.

To run a subset of unit test, use `ctest -R <test-match-pattern>` to perform tests with name matched by given pattern.

## Build with CUDA support

### Extra prerequisites

- [CUDA-Toolkit](https://developer.nvidia.com/cuda-toolkit)

To build NVIDIA GPU support for ABACUS, define `USE_CUDA` flag. You can also specify path to local installation of CUDA Toolkit by setting `CMAKE_CUDA_COMPILER` flags.

```bash
cmake -B build -DUSE_CUDA=1 -DCMAKE_CUDA_COMPILER=${path to cuda toolkit}/bin/nvcc
```

## Build math library from source

> Note: This flag is **enabled by default**. It will get better performance than the standard implementation on `gcc` and `clang`. But it **will be disabled** when using `Intel Compiler` since the math functions will get wrong results and the performance is also unexpectly poor.

To build math functions from source code, instead of using c++ standard implementation, define `USE_ABACUS_LIBM` flag. 

Currently supported math functions:
 `sin`, `cos`, `sincos`, `exp`, `cexp`

```bash
cmake -B build -DUSE_ABACUS_LIBM=1
```

## Build ABACUS with make

> Note: We suggest using CMake to configure and compile.

To compile the ABACUS program using legacy `make`, users need to specify the location of the compilers, headers and libraries in `source/Makefile.vars`:

```makefile
# This is the Makefile of ABACUS API
#======================================================================
# Users set
#======================================================================
CXX = mpiicpc
# mpiicpc:   compile intel parallel version
# icpc:      compile intel sequential version
# make: ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: nothing need to be set except LIBXC_DIR
#
# mpicxx:    compile gnu parallel version
# g++:       compile gnu sequential version
# make: FFTW_DIR, OPENBLAS_LIB_DIR, SCALAPACK_LIB_DIR, ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: FFTW_DIR, OPENBLAS_LIB_DIR must be set.

# GPU = OFF  #We do not support GPU yet
# OFF:  do not use GPU
# CUDA: use CUDA
#======================================================================



#--------------------  FOR INTEL COMPILER  ----------------------------
## ELPA_DIR          should contain an include folder and lib/libelpa.a
## CEREAL_DIR        should contain an include folder.
#----------------------------------------------------------------------

ELPA_DIR      = /usr/local/include/elpa-2021.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/elpa

CEREAL_DIR    = /usr/local/include/cereal


##-------------------  FOR GNU COMPILER  ------------------------------
## FFTW_DIR          should contain lib/libfftw3.a.
## OPENBLAS_LIB_DIR  should contain libopenblas.a. 
## SCALAPACK_LIB_DIR should contain libscalapack.a
## All three above will only be used when CXX=mpicxx or g++
## ELPA_DIR          should contain an include folder and lib/libelpa.a
## CEREAL_DIR        should contain an include folder.
##---------------------------------------------------------------------

# FFTW_DIR = /public/soft/fftw_3.3.8
# OPENBLAS_LIB_DIR   = /public/soft/openblas/lib
# SCALAPACK_LIB_DIR  = /public/soft/openblas/lib

# ELPA_DIR      = /public/soft/elpa_21.05.002
# ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002

# CEREAL_DIR    = /public/soft/cereal


##-------------------  OPTIONAL LIBS  ---------------------------------
## To use DEEPKS: set LIBTORCH_DIR and LIBNPY_DIR
## To use LIBXC:  set LIBXC_DIR which contains include and lib/libxc.a (>5.1.7)
## To use DeePMD: set DeePMD_DIR and TensorFlow_DIR
## To use LibRI:  set LIBRI_DIR and LIBCOMM_DIR
##---------------------------------------------------------------------

# LIBTORCH_DIR  = /usr/local
# LIBNPY_DIR    = /usr/local

# LIBXC_DIR    		= /public/soft/libxc

# DeePMD_DIR = ${deepmd_root}
# TensorFlow_DIR = ${tensorflow_root}

# LIBRI_DIR     = /public/software/LibRI
# LIBCOMM_DIR   = /public/software/LibComm

##---------------------------------------------------------------------
# NP = 14 # It is not supported. use make -j14 or make -j to parallelly compile
# DEBUG = OFF
# Only for developers
# ON:   use gnu compiler and check segmental defaults
# OFF:  nothing
#======================================================================
```

For example, below is a case where the Intel C++ compiler, Intel MPI and CEREAL are used, along with Intel MKL library. The file Makefile.vars can be set as
follows:

```makefile
CXX = mpiicpc #(or CXX = icpc)
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR    = /public/soft/cereal
```

When `CXX=mpiicpc`, a parallel version will be compiled. When `CXX=icpc`, a sequential version will be compiled.

Another example is where the Gnu C++ compiler, MPICH, OPENBLAS, ScaLAPACK, ELPA and CEREAL are used:

```makefile
CXX = mpicxx/g++
FFTW_DIR = /public/soft/fftw_3.3.8
OPENBLAS_LIB_DIR   = /public/soft/openblas/lib
SCALAPACK_LIB_DIR  = /public/soft/openblas/lib
ELPA_DIR      = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR    = /public/soft/cereal
```

When `CXX=mpicxx`, a parallel version will be compiled. When `CXX=g++`, a sequential version will be compiled.

Except modifying `Makefile.vars`, you can also directly use

```makefile
make CXX=mpiicpc ELPA_DIR=/public/soft/elpa_21.05.002 \
ELPA_INCLUDE_DIR=${ELPA_DIR}/include/elpa-2021.05.002 \
CEREAL_DIR=/public/soft/cereal
```

ABACUS now support full version and pw version. Use `make` or `make abacus` to compile full version which supports LCAO calculations. Use `make pw` to compile pw version which only supports pw calculations. For pw version, `make pw CXX=mpiicpc`, you do not need to provide any libs. For `make pw CXX=mpicxx`, you need provide `FFTW_DIR` and `OPENBLAS_LIB_DIR`.

Besides, libxc and deepks are optional libs to compile abacus.
They will be used when `LIBXC_DIR` is defined like

```makefile
LIBXC_DIR = /public/soft/libxc
```

or `LIBTORCH_DIR` and `LIBNPY_DIR` like

```makefile
LIBTORCH_DIR  = /usr/local
LIBNPY_DIR    = /usr/local
```

After modifying the `Makefile.vars` file, execute `make` or `make -j12` to build the program.

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

### Add DeePMD-kit Support

> Note: This part is only required if you want to load a trained DeeP Potential and run molecular dynamics with that. To train the DeeP Potential with DP-GEN, no extra prerequisite is needed and please refer to [this page](http://abacus.deepmodeling.com/en/latest/advanced/interface/dpgen.html) for ABACUS interface with DP-GEN.

To compile ABACUS with DeePMD-kit, you need to define `DeePMD_DIR` and `TensorFlow_DIR` in the file `Makefile.vars` or use

```makefile
make DeePMD_DIR=~/deepmd-kit TensorFlow_DIR=~/tensorflow
```

directly.

> `deepmd_c`/`deepmd_cc` and `tensorflow_cc` libraries would be called according to `DeePMD_DIR` and `TensorFlow_DIR`, which is showed in detail in [this page](https://github.com/deepmodeling/deepmd-kit/blob/master/doc/inference/cxx.md).

### Add LibRI and LibComm Support
To use new EXX, you need two libraries: LibRI and LibComm and need to define `LIBRI_DIR` and `LIBCOMM_DIR` in the file `Makefile.vars` or use 
```makefile
make LIBRI_DIR=/public/software/LibRI LIBCOMM_DIR=/public/software/LibComm
```
directly.
