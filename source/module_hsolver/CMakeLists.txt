list(APPEND objects
    diago_cg.cpp
    diago_david.cpp
    hsolver_pw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      hsolver_lcao.cpp
      diago_blas.cpp
      diago_elpa.cpp
      diago_lapack.cpp
  )
endif()

add_library(
    hsolver
    OBJECT
    ${objects}
)

if(ENABLE_LCAO)
  add_subdirectory(genelpa)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
