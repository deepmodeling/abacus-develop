# CMake build is an experimental feature and under active development!

Currently, only MPI version is offered.

Please build ABACUS binary with the following Bash commands:

```bash
mkdir build
cd build

cmake ../cmake
make -j 16
```

Then, the binary `ABACUS-${PROJECT_VERSION}` is generated in your working directory.

