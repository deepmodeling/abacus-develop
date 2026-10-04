abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_EXX_SYMMETRY_rotation
  LIBS base device symmetry neighbor parameter
  SOURCES symm_rotation_test.cpp test_symm_rotation.cpp ../symm_rotation.cpp ../symm_rot_out.cpp
  ../../../../source_basis/module_ao/parallel_orbitals.cpp

)