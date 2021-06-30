# Build and install ABACUS from CMake

Check the cmake version on your machine
```bash
cmake --version
```
ABACUS requires the minimum cmake version `3.18`.

Now goto the code directory `ABACUS.develop` and make a build place.
```bash
cd ABACUS.develop
mkdir build
cd build
```
You can specify the bin path of ABACUS binary to install by `CMAKE_INSTALL_PREFIX`.
```bash
cmake -DCMAKE_INSTALL_PREFIX=${ABACUS_BIN_PATH} ../cmake
```
You can provide root path of each dependent package if the package cannot be automatically found by cmake. 
Keys `ELPA_ROOT`, `FFTW3_ROOT`, `CEREAL_ROOT` and `LAPACK_ROOT` are currently available to specify. 
For example
```bash
cmake -DFFTW3_ROOT=/opt/fftw3 ../cmake
```
If the cmake has executed successfully, then
```bash
make -j 16
make install
```
If no install prefix is specified, the binary will be installed to `/usr/local/bin/ABACUS` by default.
