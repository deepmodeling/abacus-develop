abacus_disable_feature_definitions(__CUDA)

AddTest(
  TARGET dftu_pw_test
  LIBS base device parameter
  SOURCES dftu_pw_test.cpp
                ../../../source_pw/module_pwdft/dftu_tools_pw.cpp
)

AddTest(
  TARGET dftu_core_test
  LIBS base device
  SOURCES dftu_core_test.cpp
)

AddTest(
  TARGET dftu_operator_test
  LIBS base device
  SOURCES dftu_operator_test.cpp
)

if(ENABLE_LCAO AND ENABLE_MPI)
AddTest(
  TARGET dftu_lcao_test
  LIBS parameter psi base device container
  SOURCES dftu_lcao_test.cpp ../dftu_lcao_op.cpp ../dftu_fs.cpp
  ../../../source_pw/module_pwdft/dftu_base.cpp
  ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp
  ../../../source_hamilt/module_hcontainer/hcontainer.cpp
  ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp
  ../../module_operator_lcao/test/tmp_mocks.cpp
  ../../../source_hamilt/operator.cpp
)
endif()
