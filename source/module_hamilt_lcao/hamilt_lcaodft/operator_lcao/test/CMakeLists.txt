if(ENABLE_LCAO)
remove_definitions(-DUSE_NEW_TWO_CENTER)

AddTest(
  TARGET operator_overlap_test
  LIBS ${math_libs} psi base device numerical_atomic_orbitals container
  SOURCES test_overlapnew.cpp ../overlap_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp ../../../../module_hamilt_general/operator.cpp
)

AddTest(
  TARGET operator_overlap_cd_test
  LIBS ${math_libs} psi base device numerical_atomic_orbitals container
  SOURCES test_overlapnew_cd.cpp ../overlap_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp ../../../../module_hamilt_general/operator.cpp
)

AddTest(
  TARGET operator_ekinetic_test
  LIBS ${math_libs} psi base device numerical_atomic_orbitals container
  SOURCES test_ekineticnew.cpp ../ekinetic_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp ../../../../module_hamilt_general/operator.cpp
)

AddTest(
  TARGET operator_nonlocal_test
  LIBS ${math_libs} psi base device numerical_atomic_orbitals container
  SOURCES test_nonlocalnew.cpp ../nonlocal_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp ../../../../module_hamilt_general/operator.cpp
)

AddTest(
  TARGET operator_T_NL_cd_test
  LIBS ${math_libs} psi base device numerical_atomic_orbitals container 
  SOURCES test_T_NL_cd.cpp ../nonlocal_new.cpp ../ekinetic_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp ../../../../module_hamilt_general/operator.cpp
)

install(FILES parallel_operator_tests.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME operators_para_test
      COMMAND ${BASH} parallel_operator_tests.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET operator_sc_lambda_test
  LIBS ${math_libs} psi base device
  SOURCES test_sc_lambda_lcao.cpp
  ../sc_lambda_lcao.cpp
  ../../../../module_hamilt_general/operator.cpp
  ../../../module_deltaspin/spin_constrain.cpp
  ../../../module_deltaspin/cal_h_lambda.cpp
  ../../../../module_basis/module_ao/parallel_2d.cpp
  ../../../../module_basis/module_ao/parallel_orbitals.cpp
)

endif() 