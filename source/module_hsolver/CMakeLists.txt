list(APPEND objects
    diago_cg.cpp
    diago_david.cpp
    hsolver_pw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
    src/math_kernel.cpp
)

if(ENABLE_LCAO)
  list(APPEND objects
      hsolver_lcao.cpp
      diago_blas.cpp
      diago_elpa.cpp
      diago_lapack.cpp
  )
endif()

if (USE_CUDA)
  list(APPEND objects src/cuda/math_kernel.cu)
elseif(USE_ROCM)
  list(APPEND objects src/rocm/math_kernel.cu)
endif()

add_library(
    hsolver
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hsolver)
endif()

if(ENABLE_LCAO)
  add_subdirectory(genelpa)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
