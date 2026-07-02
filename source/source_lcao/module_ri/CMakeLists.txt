if (ENABLE_LIBRI)

  add_subdirectory(module_exx_symmetry)

  list(APPEND objects
        Matrix_Orbs11.cpp
        Matrix_Orbs21.cpp
        Matrix_Orbs22.cpp
        RI_2D_Comm.cpp
        Mix_DMk_2D.cpp
        Mix_Matrix.cpp
        sternheimer_abacus_dft_zero_order.cpp
        sternheimer_abacus_dft_reference.cpp
        sternheimer_abacus_fd_bindings.cpp
        sternheimer_abacus_fd_adapter.cpp
        sternheimer_dft_zero_order.cpp
        sternheimer_fd_hamiltonian.cpp
        sternheimer_fd_solver.cpp
        sternheimer_rpa.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        conv_coulomb_pot_k.cpp
        exx_abfs_ctor_orbs.cpp
        exx_abfs_io.cpp
        exx_abfs_jle.cpp
        abfs_construct_pca.cpp
        exx_opt_orb_print.cpp
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

  if(BUILD_TESTING)
    if(ENABLE_MPI)
      add_subdirectory(test)
    endif()
  endif()

  if(ENABLE_COVERAGE)
    add_coverage(ri)
  endif()
endif()
