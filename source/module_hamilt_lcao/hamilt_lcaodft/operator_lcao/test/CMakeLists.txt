if(ENABLE_LCAO)

AddTest(
  TARGET test_overlap_new
  LIBS ${math_libs} psi base device
  SOURCES test_overlapnew.cpp ../overlap_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp
)

AddTest(
  TARGET test_overlap_new_cd
  LIBS ${math_libs} psi base device
  SOURCES test_overlapnew_cd.cpp ../overlap_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp
)

AddTest(
  TARGET test_ekinetic_new
  LIBS ${math_libs} psi base device
  SOURCES test_ekineticnew.cpp ../ekinetic_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp
)

AddTest(
  TARGET test_nonlocal_new
  LIBS ${math_libs} psi base device
  SOURCES test_nonlocalnew.cpp ../nonlocal_new.cpp ../../../module_hcontainer/func_folding.cpp 
  ../../../module_hcontainer/base_matrix.cpp ../../../module_hcontainer/hcontainer.cpp ../../../module_hcontainer/atom_pair.cpp  
  ../../../../module_basis/module_ao/parallel_2d.cpp ../../../../module_basis/module_ao/parallel_orbitals.cpp 
  tmp_mocks.cpp
)

install(FILES parallel_operator_tests.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME test_operators_para
      COMMAND ${BASH} parallel_operator_tests.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

endif() 