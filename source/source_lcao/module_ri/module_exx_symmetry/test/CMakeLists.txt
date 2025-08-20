remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
AddTest(
  TARGET MODULE_RI_EXX_SYMMETRY_rotation
  LIBS base ${math_libs} device symmetry neighbor parameter
  SOURCES symmetry_rotation_test.cpp ../symmetry_rotation.cpp ../symmetry_rotation_output.cpp ../irreducible_sector.cpp ../irreducible_sector_bvk.cpp
  ../../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../../source_io/output.cpp
)