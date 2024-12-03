if(ENABLE_DEEPKS)
  list(APPEND objects
      LCAO_deepks.cpp
      deepks_fgamma.cpp
      deepks_fk.cpp
      LCAO_deepks_odelta.cpp
      LCAO_deepks_io.cpp
      LCAO_deepks_mpi.cpp
      LCAO_deepks_pdm.cpp
      LCAO_deepks_psialpha.cpp
      LCAO_deepks_torch.cpp
      LCAO_deepks_vdelta.cpp
      deepks_hmat.cpp
      LCAO_deepks_interface.cpp
      orbital_precalc.cpp
      orbital_precalc_k.cpp
      cal_gdmx.cpp
      cal_gdmx_k.cpp
      cal_gedm.cpp
      cal_gvx.cpp      
      cal_descriptor.cpp
      v_delta_precalc.cpp
      v_delta_precalc_k.cpp
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

