set(DELTA_ST_OBJECTS
  sternheimer_abacus_dft_reference.cpp
  sternheimer_abacus_dft_zero_order.cpp
  sternheimer_abacus_fd_adapter.cpp
  sternheimer_abacus_fd_bindings.cpp
  sternheimer_abacus_fd_nonlocal.cpp
  sternheimer_abacus_fd_smoke.cpp
  sternheimer_abacus_st_smoke.cpp
  sternheimer_abfs_perturbation.cpp
  sternheimer_abfs_strict2d.cpp
  sternheimer_channel_resources.cpp
  sternheimer_chi0_mpi.cpp
  sternheimer_coulomb_whitening.cpp
  sternheimer_delta.cpp
  sternheimer_dft_zero_order.cpp
  sternheimer_fd_hamiltonian.cpp
  sternheimer_fd_nonlocal_projector.cpp
  sternheimer_fd_preconditioner.cpp
  sternheimer_fd_projector_sampler.cpp
  sternheimer_fd_solver.cpp
  sternheimer_galerkin_audit.cpp
  sternheimer_galerkin_operator.cpp
  sternheimer_grid_diagnostics.cpp
  sternheimer_grid_transfer.cpp
  sternheimer_kq.cpp
  sternheimer_periodic_solver.cpp
  sternheimer_response_grid.cpp
  sternheimer_rpa.cpp
  sternheimer_runtime_options.cpp
  sternheimer_siab_memory.cpp
  sternheimer_siab_mpi.cpp
  sternheimer_siab_overlap.cpp
  sternheimer_siab_provenance.cpp
  sternheimer_siab_writer.cpp
  sternheimer_supercell_perturbation.cpp
  sternheimer_supercell_sector.cpp
  sternheimer_wavefunction_diagnostic.cpp
  sternheimer_weak_augmented.cpp
  sternheimer_weak_grid.cpp
  sternheimer_weak_preconditioner.cpp
  sternheimer_weak_q_reference.cpp
)

if(ENABLE_GREENX_MINIMAX)
  list(APPEND DELTA_ST_OBJECTS greenx_minimax_adapter.f90)
endif()

add_library(delta_st OBJECT ${DELTA_ST_OBJECTS})
target_include_directories(delta_st
  PRIVATE
    ${ABACUS_SOURCE_DIR}
    ${ABACUS_SOURCE_DIR}/source_lcao/module_ri
    ${ABACUS_SOURCE_DIR}/source_lcao/module_rpa
    ${ABACUS_SOURCE_DIR}/source_lcao/module_rpa/module_delta_st
    ${ABACUS_SOURCE_DIR}/source_lcao/module_rpa/module_librpa_io
)
target_link_libraries(delta_st PRIVATE ri rpa)

if(ENABLE_GREENX_MINIMAX)
  target_compile_definitions(delta_st PRIVATE __GREENX_MINIMAX)
  target_link_libraries(delta_st PRIVATE GreenX::MiniMax)
endif()

if(ENABLE_COVERAGE)
  add_coverage(delta_st)
endif()

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(test)
endif()
