abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES dm_mixing_test.cpp ../mix_dmk_2d.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES ri_cv_io_test.cpp
)
add_executable(MODULE_RI_ewald_tail_utils_test ewald_tail_utils_test.cpp)
install(TARGETS MODULE_RI_ewald_tail_utils_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_tail_utils_test
  COMMAND ./MODULE_RI_ewald_tail_utils_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_cutoff_test gaussian_abfs_cutoff_test.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_cutoff_test
  COMMAND ./MODULE_RI_gaussian_abfs_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_2d_integrals_test gaussian_abfs_2d_integrals_test.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_2d_integrals_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_2d_integrals_test
  COMMAND ./MODULE_RI_gaussian_abfs_2d_integrals_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_ewald_gaussian_cutoff_test ewald_gaussian_cutoff_test.cpp)
install(TARGETS MODULE_RI_ewald_gaussian_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_gaussian_cutoff_test
  COMMAND ./MODULE_RI_ewald_gaussian_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_librpa_stru_units_test librpa_stru_units_test.cpp)
install(TARGETS MODULE_RI_librpa_stru_units_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_librpa_stru_units_test
  COMMAND ./MODULE_RI_librpa_stru_units_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_librpa_bz_sampling_test librpa_bz_sampling_test.cpp)
install(TARGETS MODULE_RI_librpa_bz_sampling_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_librpa_bz_sampling_test
  COMMAND ./MODULE_RI_librpa_bz_sampling_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_librpa_2d_coulomb_head_test librpa_2d_coulomb_head_test.cpp)
install(TARGETS MODULE_RI_librpa_2d_coulomb_head_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_librpa_2d_coulomb_head_test
  COMMAND ./MODULE_RI_librpa_2d_coulomb_head_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_2d_coulomb_test
  direct_2d_coulomb_test.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_2d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_2d_coulomb_test base ${math_libs} device parameter)
install(TARGETS MODULE_RI_direct_2d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_2d_coulomb_test
  COMMAND ./MODULE_RI_direct_2d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_3d_coulomb_test
  direct_3d_coulomb_test.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_3d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_3d_coulomb_test base ${math_libs} device parameter)
install(TARGETS MODULE_RI_direct_3d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_3d_coulomb_test
  COMMAND ./MODULE_RI_direct_3d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
if(ENABLE_MPI)
  add_executable(MODULE_RI_ewald_mpi_distribution_test ewald_mpi_distribution_test.cpp)
  target_link_libraries(MODULE_RI_ewald_mpi_distribution_test MPI::MPI_CXX)
  install(TARGETS MODULE_RI_ewald_mpi_distribution_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(
    NAME MODULE_RI_ewald_mpi_distribution_test_parallel
    COMMAND mpirun -np 4 ./MODULE_RI_ewald_mpi_distribution_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

set(STERNHEIMER_RPA_TEST_LIBS ${math_libs})
if(TARGET abacus_greenx_minimax)
  list(APPEND STERNHEIMER_RPA_TEST_LIBS abacus_greenx_minimax)
endif()
AddTest(
  TARGET MODULE_RI_sternheimer_runtime_options_test
  LIBS parameter base ${math_libs} device
  SOURCES sternheimer_runtime_options_test.cpp ../sternheimer_runtime_options.cpp
)
AddTest(
  TARGET MODULE_RI_rpa_abfs_preorthogonalization_test
  LIBS parameter base ${math_libs} device
  SOURCES rpa_abfs_preorthogonalization_test.cpp
          ../rpa_abfs_preorthogonalization.cpp
          ../../../source_basis/module_ao/ORB_atomic_lm.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_rpa_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_rpa_test.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_wavefunction_diagnostic_test
  LIBS ${math_libs}
  SOURCES sternheimer_wavefunction_diagnostic_test.cpp ../sternheimer_wavefunction_diagnostic.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_supercell_perturbation_test
  LIBS ${math_libs}
  SOURCES sternheimer_supercell_perturbation_test.cpp ../sternheimer_supercell_perturbation.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_supercell_sector_test
  LIBS ${math_libs}
  SOURCES sternheimer_supercell_sector_test.cpp ../sternheimer_supercell_sector.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_hamiltonian_test
  LIBS ${math_libs}
  SOURCES sternheimer_fd_hamiltonian_test.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_preconditioner_test
  LIBS ${math_libs}
  SOURCES sternheimer_fd_preconditioner_test.cpp ../sternheimer_fd_preconditioner.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_solver_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_fd_solver_test.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_preconditioner.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_delta_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_delta_test.cpp ../sternheimer_delta.cpp ../sternheimer_kq.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_preconditioner.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_rpa.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_grid_diagnostics_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_grid_diagnostics_test.cpp ../sternheimer_grid_diagnostics.cpp ../sternheimer_delta.cpp ../sternheimer_kq.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_preconditioner.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_rpa.cpp ../../../source_base/module_external/blas_connector_matrix.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_dft_zero_order_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_dft_zero_order_test.cpp ../sternheimer_dft_zero_order.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_preconditioner.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_nonlocal_projector_test
  LIBS ${math_libs}
  SOURCES sternheimer_fd_nonlocal_projector_test.cpp ../sternheimer_fd_nonlocal_projector.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_projector_sampler_test
  LIBS parameter base ${math_libs} device container
  SOURCES sternheimer_fd_projector_sampler_test.cpp ../sternheimer_fd_projector_sampler.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abacus_fd_adapter_test
  LIBS parameter base ${math_libs} device container
  SOURCES sternheimer_abacus_fd_adapter_test.cpp ../sternheimer_abacus_fd_adapter.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_response_grid_test
  LIBS parameter base ${math_libs} device container planewave
  SOURCES sternheimer_response_grid_test.cpp ../sternheimer_response_grid.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abacus_fd_nonlocal_test
  LIBS parameter base ${math_libs} device container
  SOURCES sternheimer_abacus_fd_nonlocal_test.cpp
          ../sternheimer_abacus_fd_nonlocal.cpp
          ../sternheimer_fd_projector_sampler.cpp
          ../sternheimer_fd_nonlocal_projector.cpp
          ../sternheimer_kq.cpp
          ../../../source_cell/atom_spec.cpp
          ../../../source_cell/atom_pseudo.cpp
          ../../../source_cell/pseudo.cpp
          ../../../source_cell/test/support/mock_unitcell.cpp
          ../../../source_io/module_output/output.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abacus_fd_smoke_test
  LIBS ${math_libs}
  SOURCES sternheimer_abacus_fd_smoke_test.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abfs_perturbation_test
  LIBS parameter base ${math_libs} device container
  SOURCES sternheimer_abfs_perturbation_test.cpp ../sternheimer_abfs_perturbation.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_coulomb_tail_test
  LIBS parameter base ${math_libs} device container
  SOURCES sternheimer_coulomb_tail_test.cpp ../sternheimer_abfs_perturbation.cpp
          ../sternheimer_kq.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp
          ../conv_coulomb_pot_k.cpp ../../../source_basis/module_ao/ORB_atomic_lm.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abacus_st_smoke_test
  LIBS ${math_libs}
  SOURCES sternheimer_abacus_st_smoke_test.cpp ../sternheimer_siab_memory.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_channel_resources_test
  LIBS ${math_libs}
  SOURCES sternheimer_channel_resources_test.cpp ../sternheimer_channel_resources.cpp
)
add_executable(
  MODULE_RI_sternheimer_owned_channel_workers_test
  sternheimer_owned_channel_workers_test.cpp
  ../sternheimer_channel_resources.cpp
)
if(ENABLE_MPI)
  target_link_libraries(MODULE_RI_sternheimer_owned_channel_workers_test MPI::MPI_CXX)
endif()
install(TARGETS MODULE_RI_sternheimer_owned_channel_workers_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_sternheimer_owned_channel_workers_test
  COMMAND ./MODULE_RI_sternheimer_owned_channel_workers_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
AddTest(
  TARGET MODULE_RI_sternheimer_periodic_solver_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_periodic_solver_test.cpp
          ../sternheimer_periodic_solver.cpp
          ../sternheimer_fd_solver.cpp
          ../sternheimer_fd_hamiltonian.cpp
          ../sternheimer_fd_preconditioner.cpp
          ../sternheimer_fd_nonlocal_projector.cpp
          ../sternheimer_delta.cpp
          ../sternheimer_kq.cpp
          ../sternheimer_rpa.cpp
          ../../../source_base/module_external/blas_connector_matrix.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_siab_overlap_test
  LIBS ${math_libs}
  SOURCES test_sternheimer_siab_overlap.cpp
          ../sternheimer_siab_overlap.cpp
          ../../../source_base/module_external/blas_connector_matrix.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_siab_writer_test
  LIBS parameter ${math_libs} base device io_input
  SOURCES test_sternheimer_siab_writer.cpp ../sternheimer_siab_writer.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_siab_provenance_test
  LIBS ${math_libs}
  SOURCES test_sternheimer_siab_provenance.cpp ../sternheimer_siab_provenance.cpp
)
if(ENABLE_MPI)
  add_executable(
    MODULE_RI_sternheimer_siab_mpi_test
    test_sternheimer_siab_mpi.cpp
    ../sternheimer_chi0_mpi.cpp
    ../sternheimer_grid_diagnostics.cpp
    ../sternheimer_rpa.cpp
    ../sternheimer_siab_mpi.cpp
    ../sternheimer_siab_overlap.cpp
    ../sternheimer_siab_writer.cpp
  )
  target_link_libraries(
    MODULE_RI_sternheimer_siab_mpi_test
    parameter
    base
    device
    io_input
    ${STERNHEIMER_RPA_TEST_LIBS}
    Threads::Threads
    GTest::gtest
  )
  if(USE_OPENMP)
    target_link_libraries(MODULE_RI_sternheimer_siab_mpi_test OpenMP::OpenMP_CXX)
  endif()
  if(ENABLE_COVERAGE)
    add_coverage(MODULE_RI_sternheimer_siab_mpi_test)
  endif()
  install(TARGETS MODULE_RI_sternheimer_siab_mpi_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(
    NAME MODULE_RI_sternheimer_siab_mpi_test
    COMMAND $<TARGET_FILE:MODULE_RI_sternheimer_siab_mpi_test>
    WORKING_DIRECTORY $<TARGET_FILE_DIR:MODULE_RI_sternheimer_siab_mpi_test>
  )
  add_test(
    NAME MODULE_RI_sternheimer_siab_mpi_test_mpi2
    COMMAND mpirun -np 2 $<TARGET_FILE:MODULE_RI_sternheimer_siab_mpi_test>
    WORKING_DIRECTORY $<TARGET_FILE_DIR:MODULE_RI_sternheimer_siab_mpi_test>
  )
endif()
AddTest(
  TARGET MODULE_RI_sternheimer_siab_primitives_test
  LIBS parameter base ${math_libs} device planewave container psi
  SOURCES test_sternheimer_siab_primitives.cpp
          ../../../source_io/module_bessel/numerical_basis.cpp
          ../../../source_io/module_bessel/bessel_basis.cpp
          ../../../source_cell/parallel_kpoints.cpp
)
target_compile_options(MODULE_RI_sternheimer_siab_primitives_test PRIVATE -U__LCAO)
add_test(
  NAME MODULE_RI_sternheimer_siab_primitives_test_mpi2
  COMMAND mpirun -np 2 ./MODULE_RI_sternheimer_siab_primitives_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
