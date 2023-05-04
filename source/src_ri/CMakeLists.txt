if (ENABLE_LIBRI)
  list(APPEND objects
      exx_lip.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        conv_coulomb_pot.cpp
        conv_coulomb_pot_k.cpp
        exx_abfs-abfs_index.cpp
        exx_abfs-construct_orbs.cpp
        exx_abfs-inverse_matrix_double.cpp
        exx_abfs-io.cpp
        exx_abfs-jle.cpp
        exx_abfs-matrix_orbs11.cpp
        exx_abfs-matrix_orbs21.cpp
        exx_abfs-matrix_orbs22.cpp
        exx_abfs-pca.cpp
        exx_opt_orb-print.cpp
        exx_opt_orb.cpp
    )
  endif()
  add_library(
      ri
      OBJECT
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(ri)
  endif()
endif()