source /opt/intel/oneapi/setvars.sh -i_mpi_library_kind=debug
cmake -B build -DCMAKE_CXX_COMPILER=icpc
cmake --build build
build/ABACUS.develop/source/main
