list(APPEND objects
    diago_cg.cpp
    diago_david.cpp
    diago_bpcg.cpp
    hsolver_pw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
    diago_cusolver.cpp
    diago_cusolver.h
)

if(ENABLE_LCAO)
  list(APPEND objects
      hsolver_lcao.cpp
      diago_blas.cpp
  )
  if (USE_ELPA)
    list(APPEND objects
        diago_elpa.cpp
    )
  endif ()

  if(USE_CUSOLVER_LCAO)
      list(APPEND objects
          ./kernels/math_kernel_op.cpp
          ./kernels/dngvd_op.cpp
          ./kernels/cuda/diag_cusolver.cu

      )

    add_library(
        diag_cusolver
        OBJECT
        ${objects}
    )
    if(ENABLE_COVERAGE)
      add_coverage(diag_cusolver)
    endif()
  endif()
endif()

add_library(
    hsolver
    OBJECT
    ${objects}
)

if(ENABLE_COVERAGE)
  add_coverage(hsolver)
endif()

if(ENABLE_LCAO AND USE_ELPA)
  add_subdirectory(genelpa)
endif()

IF (BUILD_TESTING)
  add_subdirectory(test)
  add_subdirectory(kernels/test)
endif()
