
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
  TARGET test_hcontainer_output
  LIBS base ${math_libs} device
  SOURCES output_hcontainer_test.cpp tmp_mocks.cpp
    ../base_matrix.cpp
    ../hcontainer.cpp
    ../atom_pair.cpp
    ../output_hcontainer.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
    ../../../module_io/sparse_matrix.cpp  
)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET test_hcontainer_read
  LIBS base ${math_libs} device cell_info
  SOURCES read_hcontainer_test.cpp
    ../base_matrix.cpp
    ../hcontainer.cpp
    ../atom_pair.cpp
    ../output_hcontainer.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
    ../../../module_io/sparse_matrix.cpp
    ../../../module_io/csr_reader.cpp
    ../../../module_io/file_reader.cpp
    ../../../module_io/output.cpp
)

AddTest(
  TARGET test_hcontainer_ft
  LIBS base ${math_libs} device cell_info
  SOURCES test_hcontainer_ft.cpp
    ../base_matrix.cpp
    ../hcontainer.cpp
    ../atom_pair.cpp
    ../output_hcontainer.cpp
    ../../../module_basis/module_ao/parallel_2d.cpp
    ../../../module_basis/module_ao/parallel_orbitals.cpp
    ../../../module_io/sparse_matrix.cpp
    ../../../module_io/csr_reader.cpp
    ../../../module_io/file_reader.cpp
    ../../../module_io/output.cpp
)
endif() 