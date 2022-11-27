if(ENABLE_LCAO)
  if(USE_CUSOLVER_LCAO)
      list(APPEND objects
          diag_cusolver.cu
      )

    add_library(
        pdiag
        OBJECT
        ${objects}
    )
    if(ENABLE_COVERAGE)
      add_coverage(pdiag)
    endif()
  endif()
endif()
