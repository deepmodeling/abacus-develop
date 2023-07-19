remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__UT_USE_CUDA)
remove_definitions(-D__UT_USE_ROCM)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
remove_definitions(-D_OPENMP)

AddTest(
  TARGET charge_mpi_test
  LIBS ${math_libs} psi base device planewave
  SOURCES charge_mpi_test.cpp ../module_charge/charge_mpi.cpp
)
