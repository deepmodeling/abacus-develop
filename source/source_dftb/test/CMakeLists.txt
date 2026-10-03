set(dftb_native_test_sources
  test_eigensolver.cpp
  test_charge_mixing.cpp
  test_periodic_scc.cpp
  ../charge_mixing.cpp
  ../eigensolver.cpp
  ../periodic_scc.cpp
  ../sk_matrix.cpp
  ../sk_block.cpp
  ../skf_data.cpp
)

if(ENABLE_MPI)
  # Use ABACUS's shared MPI-aware GoogleTest entry point instead of making
  # direct MPI calls from the DFTB tests.
  list(APPEND dftb_native_test_sources ../../source_base/test/mpi_test_main.cpp)
endif()

AddTest(
  TARGET MODULE_DFTB_native
  LIBS parameter base device
  SOURCES ${dftb_native_test_sources}
)
