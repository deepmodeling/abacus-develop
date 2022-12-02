if(ENABLE_LCAO)
  if(USE_CUSOLVER_LCAO)
      list(APPEND objects
          diag_cusolver.cu
      )
  endif()

  add_library(
      pdiag
      OBJECT
      pdiag_double.cpp
      diag_scalapack_gvx.cpp
      ${objects}
  )

  if(ENABLE_COVERAGE)
    add_coverage(pdiag)
  endif()

#  if(BUILD_TESTING)
#    add_subdirectory(test)
#  endif()
endif()
