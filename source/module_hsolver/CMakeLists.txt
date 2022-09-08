add_library(
    hsolver
    OBJECT
    diago_cg.cpp
    diago_david.cpp
    hsolver_pw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
    hsolver_lcao.cpp
    diago_blas.cpp
    diago_elpa.cpp
    diago_lapack.cpp
)

add_subdirectory(genelpa)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
