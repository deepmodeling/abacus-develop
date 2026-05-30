remove_definitions(-D__MLALGO)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)

if(ENABLE_LCAO)
AddTest(
  TARGET MODULE_LCAO_init_dm_from_file_test
  LIBS parameter ${math_libs} base device
  SOURCES test_init_dm_from_file.cpp tmp_mocks.cpp
    ${ABACUS_SOURCE_DIR}/source_estate/module_dm/density_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_estate/module_dm/density_matrix_io.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/base_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/atom_pair.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/read_hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/func_transfer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/func_folding.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/transfer.cpp
    ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/sparse_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/csr_reader.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/file_reader.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/output.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_dm/write_dmr.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/ucell_io.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/output_hcontainer.cpp
)

AddTest(
  TARGET MODULE_LCAO_output_hcontainer_consistency_test
  LIBS parameter ${math_libs} base device
  SOURCES test_output_hcontainer_consistency.cpp tmp_mocks.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/base_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/atom_pair.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/read_hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/output_hcontainer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/func_transfer.cpp
    ${ABACUS_SOURCE_DIR}/source_lcao/module_hcontainer/transfer.cpp
    ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/sparse_matrix.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/file_reader.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/csr_reader.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_dm/write_dmr.cpp
    ${ABACUS_SOURCE_DIR}/source_io/module_output/ucell_io.cpp
)

AddTest(
  TARGET MODULE_LCAO_init_chg_hr_error_test
  LIBS parameter ${math_libs} base device
  SOURCES test_init_chg_hr_error.cpp
)
endif()
