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
        ExxLriDetail.cpp
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
