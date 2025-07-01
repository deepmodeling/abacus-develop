if (ENABLE_LIBRI)

  add_subdirectory(module_exx_symmetry)
  
  list(APPEND objects
        Matrix_Orbs11.cpp
        Matrix_Orbs21.cpp
        Matrix_Orbs22.cpp
        RI_2D_Comm.cpp
        Mix_DMk_2D.cpp
        Mix_Matrix.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        conv_coulomb_pot_k.cpp
        exx_abfs-abfs_index.cpp
        exx_abfs-construct_orbs.cpp
        exx_abfs-io.cpp
        exx_abfs-jle.cpp
        ABFs_Construct-PCA.cpp
        exx_opt_orb-print.cpp
        exx_opt_orb.cpp
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