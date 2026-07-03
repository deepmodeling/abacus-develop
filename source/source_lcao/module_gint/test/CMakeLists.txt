abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)

if(ENABLE_LCAO)

AddTest(
  TARGET MODULE_LCAO_gint_common_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_gint_common.cpp
    tmp_mocks.cpp
    ../gint_common.cpp
    ../../module_hcontainer/base_matrix.cpp
    ../../module_hcontainer/hcontainer.cpp
    ../../module_hcontainer/atom_pair.cpp
    ../../../source_basis/module_ao/parallel_orbitals.cpp
)

AddTest(
  TARGET MODULE_LCAO_gint_precision_test
  LIBS parameter ${math_libs} base device
  SOURCES test_gint_precision.cpp
    tmp_mocks.cpp
)

endif()
