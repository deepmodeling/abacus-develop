remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET dm_mixing_test
  LIBS base ${math_libs} device
  SOURCES dm_mixing_test.cpp ../Mix_DMk_2D.cpp ../Mix_Matrix.cpp
)
AddTest(
  TARGET symmetry_rotation
  LIBS base ${math_libs} device symmetry 
  SOURCES symmetry_rotation_test.cpp ../symmetry_rotation.cpp ../../module_basis/module_ao/parallel_2d.cpp
  ../../module_io/output.cpp
)