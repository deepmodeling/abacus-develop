if(ENABLE_LCAO)
abacus_disable_feature_definitions(__FFT_TWO_CENTER)

AddTest(
  TARGET MODULE_LCAO_operator_overlap_test
  LIBS parameter psi base device container
  SOURCES test_overlap.cpp ../overlap.cpp ../overlap_fs.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_hamilt/module_hcontainer/func_transfer.cpp ../../../source_hamilt/module_hcontainer/output_hcontainer.cpp ../../../source_hamilt/module_hcontainer/transfer.cpp

  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_overlap_serial_test
  LIBS parameter psi base device container
  SOURCES test_overlap_serial.cpp ../overlap.cpp ../overlap_fs.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_hamilt/module_hcontainer/func_transfer.cpp ../../../source_hamilt/module_hcontainer/output_hcontainer.cpp ../../../source_hamilt/module_hcontainer/transfer.cpp

  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_overlap_cd_test
  LIBS parameter psi base device container
  SOURCES test_overlap_cd.cpp ../overlap.cpp ../overlap_fs.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_hamilt/module_hcontainer/func_transfer.cpp ../../../source_hamilt/module_hcontainer/output_hcontainer.cpp ../../../source_hamilt/module_hcontainer/transfer.cpp

  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_ekinetic_test
  LIBS parameter psi base device container
  SOURCES test_ekinetic.cpp ../ekinetic.cpp ../ekinetic_fs.cpp ../ekinetic_dh.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_ekinetic_serial_test
  LIBS parameter psi base device container
  SOURCES test_ekinetic_serial.cpp ../ekinetic.cpp ../ekinetic_fs.cpp ../ekinetic_dh.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_nonlocal_test
  LIBS parameter psi base device container
  SOURCES test_nonlocal.cpp ../nonlocal.cpp ../nonlocal_fs.cpp ../nonlocal_dh.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_T_NL_cd_test
  LIBS parameter psi base device container
  SOURCES test_T_NL_cd.cpp ../nonlocal.cpp ../nonlocal_fs.cpp ../nonlocal_dh.cpp ../ekinetic.cpp ../ekinetic_fs.cpp ../ekinetic_dh.cpp ../operator_fs_utils.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET MODULE_LCAO_operator_dftu_test
  LIBS parameter psi base device container
  SOURCES test_dftu.cpp ../dftu_lcao.cpp ../dftu_fs.cpp ../../../source_hamilt/module_hcontainer/func_folding.cpp
  ../../../source_hamilt/module_hcontainer/base_matrix.cpp ../../../source_hamilt/module_hcontainer/hcontainer.cpp ../../../source_hamilt/module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

install(FILES parallel_operator_tests.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME MODULE_LCAO_operators_para_test
      COMMAND ${BASH} parallel_operator_tests.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

endif() 
