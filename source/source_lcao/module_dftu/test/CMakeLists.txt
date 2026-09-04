abacus_disable_feature_definitions(__CUDA)

AddTest(
  TARGET dftu_pw_test
  LIBS base device parameter
  SOURCES dftu_pw_test.cpp
                ../../../source_pw/module_pwdft/dftu_base_tools.cpp
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
  SOURCES dftu_lcao_test.cpp ../dftu_nao_op.cpp ../dftu_nao_fs_r.cpp ../dftu_nao_for_r.cpp ../dftu_nao_str_r.cpp
  ../../../source_pw/module_pwdft/dftu_base.cpp
  ../../../source_pw/module_pwdft/dftu_base_io.cpp
  ../../../source_pw/module_pwdft/yukawa_screening.cpp
  ../../../source_estate/occ_matrix.cpp
  ../../../source_estate/occ_mixer.cpp
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
