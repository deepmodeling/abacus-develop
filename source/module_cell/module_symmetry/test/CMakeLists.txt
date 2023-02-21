remove_definitions(-D__LCAO)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET symmetry_symmetry
  LIBS base ${math_libs} device
  SOURCES symmetry_test.cpp ../symmetry.cpp ../symmetry_basic.cpp ../symm_other.cpp 
)
