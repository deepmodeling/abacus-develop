list(APPEND objects
    diag_const_nums.cpp
    diago_cg.cpp
    diago_david.cpp
    diago_dav_subspace.cpp
    diago_bpcg.cpp
    hsolver_pw.cpp
    hsolver_lcaopw.cpp
    hsolver_pw_sdft.cpp
    diago_iter_assist.cpp
    hsolver.cpp
    diago_pxxxgvx.cpp
    diag_hs_para.cpp

)

if(ENABLE_LCAO)
  if (ENABLE_MPI)
    list(APPEND objects
        hsolver_lcao.cpp
        diago_scalapack.cpp
        parallel_k2d.cpp
    )
  else ()
    list(APPEND objects
        hsolver_lcao.cpp
        diago_lapack.cpp
    )
  endif ()
  if (USE_ELPA)
    list(APPEND objects
        diago_elpa.cpp
        diago_elpa_native.cpp
    )
  endif ()

  if(USE_CUDA)
      list(APPEND objects
          ./kernels/dngvd_op.cpp
          ./kernels/cuda/diag_cusolver.cu
          diago_cusolver.cpp
          diago_cusolver.h
      )
      if(ENABLE_CUSOLVERMP)
        list(APPEND objects
            ./kernels/cuda/diag_cusolvermp.cu
            diago_cusolvermp.cpp
            diago_cusolvermp.h
        )
      endif()
    add_library(
        diag_cusolver
        OBJECT
        ${objects}
    )
    if(ENABLE_COVERAGE)
      add_coverage(diag_cusolver)
    endif()
  endif()

  if(ENABLE_PEXSI)
  list(APPEND objects
      diago_pexsi.cpp
    )
    add_subdirectory(module_pexsi)
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
  if(ENABLE_MPI)
    add_subdirectory(kernels/test)
  endif()
  message(STATUS "Building tests")
endif()
