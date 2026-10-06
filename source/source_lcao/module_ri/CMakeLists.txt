if (ENABLE_LIBRI)

  add_subdirectory(module_exx_symmetry)

  list(APPEND objects
        matrix_orbs11.cpp
        matrix_orbs21.cpp
        matrix_orbs22.cpp
        ri_2d_comm.cpp
        mix_dmk_2d.cpp
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
        rpa_lri.cpp
        rpa_lri_detail.cpp
        rpa_lri_mpi_layout.cpp
        rpa_lri_wfc_mpi.cpp
        rpa_lri_coulomb.cpp
        rpa_lri_coulomb_io.cpp
        rpa_lri_overlap.cpp
        rpa_lri_overlap_io.cpp
        rpa_lri_output.cpp
        rpa_lri_coeff_io.cpp
        rpa_lri_velocity.cpp
        rpa_lri_basis.cpp
        rpa_lri_eigenvector_v1.cpp
        rpa_lri_eigenvector.cpp
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
