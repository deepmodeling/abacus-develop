if (ENABLE_LIBRI)
  list(APPEND objects
      exx_lip.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        conv_coulomb_pot_k.cpp
        exx_abfs-abfs_index.cpp
        exx_abfs-construct_orbs.cpp
        exx_abfs-io.cpp
        exx_abfs-jle.cpp
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