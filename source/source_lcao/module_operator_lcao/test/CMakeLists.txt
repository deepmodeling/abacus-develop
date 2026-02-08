if(ENABLE_LCAO)
remove_definitions(-DUSE_NEW_TWO_CENTER)

AddTest(
  TARGET operator_overlap_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_overlap.cpp ../overlap.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../module_hcontainer/func_transfer.cpp ../../module_hcontainer/output_hcontainer.cpp ../../module_hcontainer/transfer.cpp
  ../../../source_io/sparse_matrix.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET operator_overlap_serial_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_overlap_serial.cpp ../overlap.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../module_hcontainer/func_transfer.cpp ../../module_hcontainer/output_hcontainer.cpp ../../module_hcontainer/transfer.cpp
  ../../../source_io/sparse_matrix.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET operator_overlap_cd_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_overlap_cd.cpp ../overlap.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../module_hcontainer/func_transfer.cpp ../../module_hcontainer/output_hcontainer.cpp ../../module_hcontainer/transfer.cpp
  ../../../source_io/sparse_matrix.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
  ../../module_rt/td_folding.cpp
)

AddTest(
  TARGET operator_ekinetic_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_ekinetic.cpp ../ekinetic.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET operator_ekinetic_serial_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_ekinetic_serial.cpp ../ekinetic.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET operator_nonlocal_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_nonlocal.cpp ../nonlocal.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET operator_T_NL_cd_test
  LIBS parameter ${math_libs} psi base device container
  SOURCES test_T_NL_cd.cpp ../nonlocal.cpp ../ekinetic.cpp ../operator_force_stress_utils.cpp ../../module_hcontainer/func_folding.cpp
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp
  ../../../source_basis/module_ao/parallel_orbitals.cpp
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

AddTest(
  TARGET operator_dftu_test
  LIBS parameter ${math_libs} psi base device container 
  SOURCES test_dftu.cpp ../dftu_lcao.cpp ../../module_hcontainer/func_folding.cpp 
  ../../module_hcontainer/base_matrix.cpp ../../module_hcontainer/hcontainer.cpp ../../module_hcontainer/atom_pair.cpp  
  ../../../source_basis/module_ao/parallel_orbitals.cpp 
  ../../../source_basis/module_ao/ORB_atomic_lm.cpp
  tmp_mocks.cpp ../../../source_hamilt/operator.cpp
)

install(FILES parallel_operator_tests.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME operators_para_test
      COMMAND ${BASH} parallel_operator_tests.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

endif() 
