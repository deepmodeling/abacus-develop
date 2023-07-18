if(ENABLE_LCAO)
  add_library(
      orb
      OBJECT
      ORB_atomic.cpp
      ORB_atomic_lm.cpp
      ORB_control.cpp
      ORB_gaunt_table.cpp
      ORB_gen_tables.cpp
      ORB_nonlocal.cpp
      ORB_nonlocal_lm.cpp
      ORB_read.cpp
      ORB_table_alpha.cpp
      ORB_table_beta.cpp
      ORB_table_phi.cpp
      parallel_2d.cpp
      parallel_orbitals.cpp
  )

  if(ENABLE_COVERAGE)
    add_coverage(orb)
  endif()

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
