if(ENABLE_LCAO)
  add_library(
      orb
      OBJECT
      orb_atomic.cpp
      orb_atomic_lm.cpp
      orb_gaunt_table.cpp
      orb_nonlocal.cpp
      orb_nonlocal_lm.cpp
      orb_read.cpp
      parallel_orbitals.cpp
      elem_basis_idx_orb.cpp
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
