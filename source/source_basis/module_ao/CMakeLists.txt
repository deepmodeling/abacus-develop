if(ENABLE_LCAO)
  add_library(
      orb
      OBJECT
      ORB_atomic.cpp
      ORB_atomic_lm.cpp
      ORB_gaunt_table.cpp
      ORB_nonlocal.cpp
      ORB_nonlocal_lm.cpp
      ORB_read.cpp
      parallel_orbitals.cpp
  )

  if(ENABLE_COVERAGE)
    add_coverage(orb)
  endif()

  if(BUILD_TESTING)
    if(ENABLE_MPI)
      add_subdirectory(test)
    endif()
  endif()
endif()
