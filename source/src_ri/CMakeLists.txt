if (ENABLE_LIBRI)
  list(APPEND objects
      exx_lip.cpp
  )

  if(ENABLE_LCAO)
    list(APPEND objects
        abfs.cpp
        conv_coulomb_pot.cpp
        conv_coulomb_pot_k.cpp
        exx_abfs-abfs_index.cpp
        exx_abfs-construct_orbs.cpp
        exx_abfs-dm.cpp
        exx_abfs-inverse_matrix_double.cpp
        exx_abfs-io.cpp
        exx_abfs-jle.cpp
        exx_abfs-matrix_lcaoslcaos_lcaoslcaos.cpp
        exx_abfs-matrix_orbs11.cpp
        exx_abfs-matrix_orbs21.cpp
        exx_abfs-matrix_orbs22.cpp
        # exx_abfs-parallel-communicate-allreduce.cpp
        # exx_abfs-parallel-communicate-dm-allreduce.cpp
        # exx_abfs-parallel-communicate-dm.cpp
        # exx_abfs-parallel-communicate-dm2.cpp
        exx_abfs-parallel-communicate-dm3-allreduce.cpp
        exx_abfs-parallel-communicate-dm3.cpp
        exx_abfs-parallel-communicate-function.cpp
        # exx_abfs-parallel-communicate-hexx-allreduce.cpp
        exx_abfs-parallel-communicate-hexx-allreduce2.cpp
        exx_abfs-parallel-communicate-hexx.cpp
        exx_abfs-parallel-distribute-htime.cpp
        exx_abfs-parallel-distribute-kmeans.cpp
        exx_abfs-parallel-distribute-order.cpp
        exx_abfs-pca.cpp
        exx_abfs-screen-cauchy.cpp
        exx_abfs-screen-schwarz.cpp
        exx_abfs-util.cpp
        exx_abfs.cpp
        exx_lcao.cpp
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