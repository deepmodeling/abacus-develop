
if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_LCAO_hcontainer_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_hcontainer.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_LCAO_hcontainer_complex_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_hcontainer_complex.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_LCAO_hcontainer_cost_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_hcontainer_time.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_LCAO_hcontainer_folding_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_func_folding.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../func_folding.cpp ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET MODULE_LCAO_hcontainer_transfer_test
  LIBS parameter ${math_libs} psi base device
  SOURCES test_transfer.cpp ../func_transfer.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../transfer.cpp ../../../source_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

install(FILES parallel_hcontainer_tests.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME MODULE_LCAO_hcontainer_para_test
      COMMAND ${BASH} parallel_hcontainer_tests.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

AddTest(
  TARGET MODULE_LCAO_hcontainer_output_test
  LIBS parameter base ${math_libs} device
  SOURCES test_hcontainer_output.cpp
    tmp_mocks.cpp
    ../output_hcontainer.cpp
    ../base_matrix.cpp
    ../hcontainer.cpp
    ../atom_pair.cpp
    ../../../source_basis/module_ao/parallel_orbitals.cpp
    ../../../source_io/module_output/sparse_matrix.cpp
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


endif() 
