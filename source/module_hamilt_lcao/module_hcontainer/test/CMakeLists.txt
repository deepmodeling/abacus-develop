
if(ENABLE_LCAO)
AddTest(
  TARGET test_hcontainer
  LIBS ${math_libs} psi base device
  SOURCES test_hcontainer.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../module_basis/module_ao/parallel_2d.cpp ../../../module_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET test_hcontainer_complex
  LIBS ${math_libs} psi base device
  SOURCES test_hcontainer_complex.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../module_basis/module_ao/parallel_2d.cpp ../../../module_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET test_hcontainer_cost
  LIBS ${math_libs} psi base device
  SOURCES test_hcontainer_time.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp  
  ../../../module_basis/module_ao/parallel_2d.cpp ../../../module_basis/module_ao/parallel_orbitals.cpp tmp_mocks.cpp
)

AddTest(
  TARGET test_output_hcontainer
  LIBS base ${math_libs} device
  SOURCES output_hcontainer_test.cpp ../output_hcontainer.cpp ../../../module_basis/module_ao/parallel_2d.cpp ../../../module_basis/module_ao/parallel_orbitals.cpp ../../../module_io/sparse_matrix.cpp ../base_matrix.cpp ../hcontainer.cpp ../atom_pair.cpp tmp_mocks.cpp
)
endif() 