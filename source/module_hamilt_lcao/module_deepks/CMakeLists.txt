if(ENABLE_DEEPKS)
  list(APPEND objects
      LCAO_deepks.cpp
      deepks_force.cpp
      deepks_orbital.cpp
      LCAO_deepks_io.cpp
      LCAO_deepks_mpi.cpp
      LCAO_deepks_pdm.cpp
      LCAO_deepks_phialpha.cpp
      LCAO_deepks_torch.cpp
      LCAO_deepks_vdelta.cpp
      deepks_hmat.cpp
      LCAO_deepks_interface.cpp
      deepks_orbpre.cpp
      cal_gdmx.cpp
      cal_gedm.cpp
      cal_gvx.cpp      
      cal_descriptor.cpp
      v_delta_precalc.cpp
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
#  if(BUILD_TESTING)
#      add_subdirectory(test)
#  endif()
endif()

