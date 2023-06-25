if(ENABLE_DEEPKS)
  list(APPEND objects
      LCAO_deepks.cpp
      LCAO_deepks_fdelta.cpp
      LCAO_deepks_odelta.cpp
      LCAO_deepks_io.cpp
      LCAO_deepks_mpi.cpp
      LCAO_deepks_pdm.cpp
      LCAO_deepks_psialpha.cpp
      LCAO_deepks_torch.cpp
      LCAO_deepks_vdelta.cpp
      LCAO_deepks_interface.cpp
  )

  add_library(
      deepks
      OBJECT
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(deepks)
  endif()

  if(BUILD_TESTING)
      add_subdirectory(test)
  endif()
endif()

