if(ENABLE_LCAO)
  add_subdirectory(MRRR)
  if(USE_CUSOLVER_LCAO)
      list(APPEND objects
          diag_cusolver.cu
      )
  endif()

  add_library(
      pdiag
      OBJECT
      pdgseps.cpp
      pdiag_common.cpp
      pdiag_double.cpp
      pdst2g.cpp
      pdstebz.cpp
      pdsteiz.cpp
      pdsyg2st.cpp
      pdsytrd.cpp
      pdt2s.cpp
      pdtrsm.cpp
      pzgseps.cpp
      pzheg2st.cpp
      pzhetrd.cpp
      pzhtrsm.cpp
      pzst2g.cpp
      pzsteiz.cpp
      pzt2s.cpp
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
