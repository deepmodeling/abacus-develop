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

set(STERNHEIMER_RPA_TEST_LIBS ${math_libs})
if(TARGET abacus_greenx_minimax)
  list(APPEND STERNHEIMER_RPA_TEST_LIBS abacus_greenx_minimax)
endif()
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
  TARGET MODULE_RI_sternheimer_fd_solver_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_fd_solver_test.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_delta_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_delta_test.cpp ../sternheimer_delta.cpp ../sternheimer_kq.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_dft_zero_order_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_dft_zero_order_test.cpp ../sternheimer_dft_zero_order.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_fd_nonlocal_projector.cpp ../sternheimer_rpa.cpp
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
  TARGET MODULE_RI_sternheimer_abacus_st_smoke_test
  LIBS ${math_libs}
  SOURCES sternheimer_abacus_st_smoke_test.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_periodic_solver_test
  LIBS ${STERNHEIMER_RPA_TEST_LIBS}
  SOURCES sternheimer_periodic_solver_test.cpp
          ../sternheimer_periodic_solver.cpp
          ../sternheimer_fd_solver.cpp
          ../sternheimer_fd_hamiltonian.cpp
          ../sternheimer_fd_nonlocal_projector.cpp
          ../sternheimer_delta.cpp
          ../sternheimer_kq.cpp
          ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_kq_test
  LIBS ${math_libs}
  SOURCES sternheimer_kq_test.cpp ../sternheimer_kq.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_channel_resources_test
  LIBS ${math_libs}
  SOURCES sternheimer_channel_resources_test.cpp ../sternheimer_channel_resources.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_siab_overlap_test
  LIBS ${math_libs}
  SOURCES test_sternheimer_siab_overlap.cpp ../sternheimer_siab_overlap.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_siab_writer_test
  LIBS parameter ${math_libs} base device io_input
  SOURCES test_sternheimer_siab_writer.cpp ../sternheimer_siab_writer.cpp
)
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
