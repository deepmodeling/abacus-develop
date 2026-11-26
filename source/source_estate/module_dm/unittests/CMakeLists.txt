abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

abacus_copy_test_files(support)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_ESTATE_edm_tddft
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter base container device
    SOURCES test_edm_tddft.cpp ../edm_tddft.cpp ../edm_tddft_lapack.cpp
            ../../../source_base/test/mpi_test_main.cpp)
  add_test(
    NAME MODULE_ESTATE_edm_tddft_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_ESTATE_edm_tddft> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_ESTATE_edm_tddft MODULE_ESTATE_edm_tddft_np4
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
  set_tests_properties(MODULE_ESTATE_edm_tddft_np4 PROPERTIES PROCESSORS 4)

  addtest(
    TARGET MODULE_ESTATE_dm_constructor_test
    LIBS parameter base device symmetry
    SOURCES test_dm_constructor.cpp
            ../density_matrix.cpp
            ../dmr_gamma.cpp
            ../dmr_init.cpp
            ../dm_setter.cpp
            ../dm_getter.cpp
            tmp_mocks.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
            ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp)

  addtest(
    TARGET MODULE_ESTATE_dm_init_test
    LIBS parameter base device symmetry
    SOURCES test_dm_r_init.cpp
            ../density_matrix.cpp
            ../dmr_gamma.cpp
            ../dmr_init.cpp
            ../dm_setter.cpp
            ../dm_getter.cpp
            tmp_mocks.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
            ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
            ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp)

endif()

# The calculation fixture initializes both serial and distributed orbital layouts.
addtest(
  TARGET MODULE_ESTATE_dm_cal_DMR_test
  LIBS parameter base device symmetry
  SOURCES test_cal_dm_r.cpp
          ../density_matrix.cpp
          ../dmr_gamma.cpp
          ../dmr_init.cpp
          ../dm_setter.cpp
          ../dm_getter.cpp
          ../dm_tools.cpp
          ../dmr_k.cpp
          ../dmr_td.cpp
          ../dmr_full.cpp
          tmp_mocks.cpp
          ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
          ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
          ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
          ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp
          ${ABACUS_SOURCE_DIR}/source_cell/klist.cpp
          ${ABACUS_SOURCE_DIR}/source_cell/klist_io.cpp
          ${ABACUS_SOURCE_DIR}/source_cell/parallel_kpoints.cpp
          ${ABACUS_SOURCE_DIR}/source_cell/reciprocal_grid.cpp)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_ESTATE_dm_soc_magnetization_roundtrip_test
    LIBS parameter base device
    SOURCES test_soc_magnetization_roundtrip.cpp
            ../density_matrix.cpp
            ../dmr_gamma.cpp
            ../dmr_init.cpp
            ../dm_setter.cpp
            ../dm_getter.cpp
            ../dm_tools.cpp
            ../dmr_k.cpp
            ../dmr_td.cpp
            ../dmr_full.cpp
            tmp_mocks.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/base_matrix.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/hcontainer.cpp
            ${ABACUS_SOURCE_DIR}/source_hamilt/module_hcontainer/atom_pair.cpp
            ${ABACUS_SOURCE_DIR}/source_basis/module_ao/parallel_orbitals.cpp)
endif()
