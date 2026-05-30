remove_definitions(-D__MPI)
remove_definitions(-D__CUDA)
remove_definitions(-D__UT_USE_CUDA)
remove_definitions(-D__UT_USE_ROCM)
remove_definitions(-D__ROCM)

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
