# CMake build is an experimental feature and under active development!

Currently, only MPI version is offered.

Please build ABACUS binary with the following Bash commands:

```bash
mkdir build
cd build

cmake ../cmake

# Or specifying installation paths of dependencies:
# cmake ../cmake \
#   -DCEREAL_DIR=${CEREAL_DIR} \
#   -DELPA_DIR=${FFTW_DIR} \
#   -DFFTW_DIR=${ELPA_DIR} \
#   -DSCALAPACK_DIR=${SCALAPACK_DIR}

make -j 16
```

Then, the binary `ABACUS-${PROJECT_VERSION}` is generated in your working directory.

