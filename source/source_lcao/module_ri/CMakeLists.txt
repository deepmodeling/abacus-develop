if (ENABLE_LIBRI)

  add_subdirectory(module_exx_symmetry)
  if(ENABLE_GREENX_MINIMAX)
    add_subdirectory(greenx_minimax)
  endif()
  
  list(APPEND objects
        Matrix_Orbs11.cpp
        Matrix_Orbs21.cpp
        Matrix_Orbs22.cpp
        RI_2D_Comm.cpp
        Mix_DMk_2D.cpp
        Mix_Matrix.cpp
        rpa_abfs_preorthogonalization.cpp
        direct_2d_coulomb.cpp
        sternheimer_abacus_dft_zero_order.cpp
        sternheimer_abacus_dft_reference.cpp
        sternheimer_abfs_perturbation.cpp
        sternheimer_abacus_fd_bindings.cpp
        sternheimer_abacus_fd_adapter.cpp
        sternheimer_abacus_fd_nonlocal.cpp
        sternheimer_abacus_fd_smoke.cpp
        sternheimer_channel_resources.cpp
        sternheimer_chi0_mpi.cpp
        sternheimer_coulomb_whitening.cpp
        sternheimer_dft_zero_order.cpp
        sternheimer_delta.cpp
        sternheimer_grid_diagnostics.cpp
        sternheimer_fd_hamiltonian.cpp
        sternheimer_fd_preconditioner.cpp
        sternheimer_fd_nonlocal_projector.cpp
        sternheimer_fd_projector_sampler.cpp
        sternheimer_fd_solver.cpp
        sternheimer_kq.cpp
        sternheimer_periodic_solver.cpp
        sternheimer_rpa.cpp
        sternheimer_response_grid.cpp
        sternheimer_runtime_options.cpp
        sternheimer_siab_mpi.cpp
        sternheimer_siab_memory.cpp
        sternheimer_siab_overlap.cpp
        sternheimer_siab_provenance.cpp
        sternheimer_siab_writer.cpp
        sternheimer_supercell_perturbation.cpp
        sternheimer_supercell_sector.cpp
        sternheimer_wavefunction_diagnostic.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        conv_coulomb_pot_k.cpp
        sternheimer_abacus_st_smoke.cpp
        exx_abfs-construct_orbs.cpp
        exx_abfs-io.cpp
        exx_abfs-jle.cpp
        ABFs_Construct-PCA.cpp
        exx_opt_orb-print.cpp
        exx_opt_orb.cpp
        gaussian_abfs.cpp
        singular_value.cpp
        exx_lri_detail.cpp
    )
  endif()
  add_library(
      ri
      OBJECT
      ${objects}
  )

  if(TARGET abacus_greenx_minimax)
    target_link_libraries(ri PUBLIC abacus_greenx_minimax)
  endif()
  
  if(BUILD_TESTING)
    if(ENABLE_MPI)
      add_subdirectory(test)
    endif()
  endif()

  if(ENABLE_COVERAGE)
    add_coverage(ri)
  endif()
endif()
