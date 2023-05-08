remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)

AddTest(
  TARGET Elecstate_Op_UTs
  LIBS ${math_libs} psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
)

AddTest(
  TARGET elecstate_magnetism
  LIBS ${math_libs} base device
  SOURCES elecstate_magnetism_test.cpp ../magnetism.cpp
)

AddTest(
  TARGET elecstate_energy
  LIBS ${math_libs} base device
  SOURCES elecstate_energy_test.cpp ../energy.cpp
)
