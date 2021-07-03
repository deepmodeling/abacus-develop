cd ABACUS.develop
cmake -B build cmake
cmake --build build -j4
cmake --install build
cd build
ctest --verbose
