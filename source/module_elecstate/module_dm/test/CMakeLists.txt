remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__MPI)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET dm_io_test
  LIBS base ${math_libs} base device cell_info
  SOURCES test_dm_io.cpp ../density_matrix.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/output_hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/module_basis/module_ao/parallel_2d.cpp
  ${ABACUS_SOURCE_DIR}/module_basis/module_ao/parallel_orbitals.cpp
  ${ABACUS_SOURCE_DIR}/module_io/sparse_matrix.cpp
  ${ABACUS_SOURCE_DIR}/module_io/csr_reader.cpp
  ${ABACUS_SOURCE_DIR}/module_io/file_reader.cpp
  ${ABACUS_SOURCE_DIR}/module_io/output.cpp
)

AddTest(
  TARGET dm_constructor_test
  LIBS base ${math_libs} base device cell_info
  SOURCES test_dm_constructor.cpp ../density_matrix.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/module_hamilt_lcao/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/module_basis/module_ao/parallel_2d.cpp
  ${ABACUS_SOURCE_DIR}/module_basis/module_ao/parallel_orbitals.cpp
  ${ABACUS_SOURCE_DIR}/module_io/sparse_matrix.cpp
  ${ABACUS_SOURCE_DIR}/module_io/csr_reader.cpp
  ${ABACUS_SOURCE_DIR}/module_io/file_reader.cpp
  ${ABACUS_SOURCE_DIR}/module_io/output.cpp
)
