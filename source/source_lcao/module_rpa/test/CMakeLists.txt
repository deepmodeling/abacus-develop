abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

AddTest(
  TARGET MODULE_RPA_rpa_abfs_preorthogonalization_test
  LIBS parameter base ${math_libs} device
  SOURCES rpa_abfs_preorthogonalization_test.cpp
          ../rpa_abfs_preorthogonalization.cpp
          ../../../source_basis/module_ao/ORB_atomic_lm.cpp
)
