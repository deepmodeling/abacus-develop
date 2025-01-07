if(ENABLE_DEEPKS)
  list(APPEND objects
      LCAO_deepks.cpp
      deepks_descriptor.cpp
      deepks_force.cpp
      deepks_fpre.cpp
      deepks_spre.cpp
      deepks_orbital.cpp
      deepks_orbpre.cpp
      deepks_vdpre.cpp
      deepks_hmat.cpp
      LCAO_deepks_io.cpp
      LCAO_deepks_pdm.cpp
      LCAO_deepks_phialpha.cpp
      LCAO_deepks_torch.cpp
      LCAO_deepks_vdelta.cpp
      LCAO_deepks_interface.cpp
      cal_gedm.cpp    
  )

  add_library(
      deepks
      OBJECT
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(deepks)
  endif()

# I will rewrite the test later, the current test rely on too many modules
  if(BUILD_TESTING)
      add_subdirectory(test)
  endif()
endif()

