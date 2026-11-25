abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)

if(ENABLE_LCAO)

  addtest(
    TARGET MODULE_LCAO_gint_common_test
    LIBS parameter psi base device
    SOURCES test_gint_common.cpp
            tmp_mocks.cpp
            ../gint_common.cpp
            ../../../source_hamilt/module_hcontainer/base_matrix.cpp
            ../../../source_hamilt/module_hcontainer/hcontainer.cpp
            ../../../source_hamilt/module_hcontainer/atom_pair.cpp
            ../../../source_basis/module_ao/parallel_orbitals.cpp)

  addtest(
    TARGET MODULE_LCAO_gint_precision_test
    LIBS parameter base device
    SOURCES test_gint_precision.cpp tmp_mocks.cpp)

  addtest(
    TARGET MODULE_LCAO_gint_prec_ctrl_test
    LIBS parameter base device
    SOURCES test_gint_prec_ctrl.cpp ../gint_prec_ctrl.cpp)

endif()
