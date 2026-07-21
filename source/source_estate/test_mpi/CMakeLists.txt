abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

AddTest(
  TARGET MODULE_ESTATE_charge_mpi_test
  LIBS parameter psi base device planewave
  SOURCES charge_mpi_test.cpp ../module_charge/charge_mpi.cpp
)

add_test(NAME MODULE_ESTATE_charge_mpi_test_4np
      COMMAND mpirun -np 4 ./MODULE_ESTATE_charge_mpi_test;
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
