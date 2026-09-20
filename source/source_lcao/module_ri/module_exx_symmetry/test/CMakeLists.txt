abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_EXX_SYMMETRY_rotation
  LIBS base device symmetry neighbor parameter
  SOURCES symm_rotation_test.cpp test_symm_rotation.cpp ../symm_rotation.cpp ../symm_rot_out.cpp
  ../../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../../source_io/module_output/output.cpp
)
AddTest(
  TARGET MODULE_RI_EXX_SYMMETRY_abf_rotation
  LIBS base ${math_libs} device symmetry neighbor parameter
  SOURCES abf_rotation_test.cpp ../symmetry_rotation.cpp ../symmetry_rotation_output.cpp ../irreducible_sector.cpp ../irreducible_sector_bvk.cpp
  ../../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../../source_io/module_output/output.cpp
)
AddTest(
  TARGET MODULE_RI_EXX_SYMMETRY_nspin4_restore
  LIBS base ${math_libs} device symmetry neighbor parameter
  SOURCES nspin4_restore_test.cpp ../symmetry_rotation.cpp ../symmetry_rotation_output.cpp ../irreducible_sector.cpp ../irreducible_sector_bvk.cpp
  ../../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../../source_io/module_output/output.cpp
)
