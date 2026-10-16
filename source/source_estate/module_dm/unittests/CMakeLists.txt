abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_ESTATE_dm_constructor_test
  LIBS parameter  base device symmetry
  SOURCES test_dm_constructor.cpp ../density_matrix.cpp ../dmr_gamma.cpp ../dmr_init.cpp ../dm_setter.cpp ../dm_getter.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_dm_init_test
  LIBS parameter  base device symmetry
  SOURCES test_dm_r_init.cpp ../density_matrix.cpp ../dmr_gamma.cpp ../dmr_init.cpp ../dm_setter.cpp ../dm_getter.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_dm_cal_DMR_test
  LIBS parameter  base device symmetry
  SOURCES test_cal_dm_r.cpp ../density_matrix.cpp ../dmr_gamma.cpp ../dmr_init.cpp ../dm_setter.cpp ../dm_getter.cpp ../dm_tools.cpp ../dmr_k.cpp ../dmr_td.cpp ../dmr_full.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_dm_soc_magnetization_roundtrip_test
  LIBS parameter  base device
  SOURCES test_soc_magnetization_roundtrip.cpp ../density_matrix.cpp ../dmr_gamma.cpp ../dmr_init.cpp ../dm_setter.cpp ../dm_getter.cpp ../dm_tools.cpp ../dmr_k.cpp ../dmr_td.cpp ../dmr_full.cpp tmp_mocks.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
  ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
  ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
)
