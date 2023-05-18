remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
remove_definitions(-D_OPENMP)

AddTest(
  TARGET Elecstate_Op_UTs
  LIBS ${math_libs} psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
)

AddTest(
  TARGET elecstate_occupy
  LIBS ${math_libs} base device
  SOURCES elecstate_occupy_test.cpp ../occupy.cpp
)

AddTest(
  TARGET elecstate_magnetism
  LIBS ${math_libs} base device
  SOURCES elecstate_magnetism_test.cpp ../magnetism.cpp
)

AddTest(
  TARGET elecstate_fp_energy
  LIBS ${math_libs} base device
  SOURCES elecstate_fp_energy_test.cpp ../fp_energy.cpp
)

AddTest(
  TARGET elecstate_print
  LIBS ${math_libs} base device
  SOURCES elecstate_print_test.cpp ../elecstate_print.cpp
)

AddTest(
  TARGET elecstate_base
  LIBS ${math_libs} base device
  SOURCES elecstate_base_test.cpp ../elecstate.cpp ../occupy.cpp ../../module_psi/psi.cpp
)
