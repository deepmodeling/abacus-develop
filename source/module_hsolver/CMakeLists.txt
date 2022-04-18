add_library(
    hsolver
    OBJECT
    diago_cg.cpp
    diago_david.cpp
    hsolver_pw.cpp
    diago_iter_assist.cpp
    hsolver_lcao.cpp
    diago_blas.cpp
    diago_elpa.cpp
)

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
