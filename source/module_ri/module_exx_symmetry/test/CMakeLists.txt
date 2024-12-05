remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET symmetry_rotation
  LIBS base ${math_libs} device symmetry neighbor parameter
  SOURCES symmetry_rotation_test.cpp ../symmetry_rotation.cpp ../irreducible_sector.cpp ../irreducible_sector_bvk.cpp
  ../../../module_basis/module_ao/parallel_orbitals.cpp
  ../../../module_io/output.cpp
)