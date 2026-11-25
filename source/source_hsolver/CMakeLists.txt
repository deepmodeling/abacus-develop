list(
  APPEND
  objects
  diag_const_nums.cpp
  diago_cg.cpp
  diago_david.cpp
  diago_dav_subspace.cpp
  diago_bpcg.cpp
  para_lin_tf.cpp
  hsolver_pw.cpp
  hsolver_lcaopw.cpp
  diago_iter_assist.cpp
  hsolver.cpp
  diago_pxxxgvx.cpp
  diag_hs_para.cpp
  diago_params.cpp)

if(ENABLE_LCAO)
  if(ENABLE_MPI)
    list(APPEND objects hsolver_lcao.cpp diago_scalapack.cpp parallel_k2d.cpp
         diago_lapack.cpp)
  else()
    list(APPEND objects hsolver_lcao.cpp diago_lapack.cpp)
  endif()
  if(ENABLE_ELPA)
    list(APPEND objects diago_elpa.cpp diago_elpa_native.cpp)
  endif()

  # DSP build without C_PACK: provide the zhegvx_ fallback for
  # cpack_zhegvx_driver (used by diago_lapack.cpp under __DSP).
  if(USE_DSP AND NOT MT_CPACK_DIR)
    list(APPEND objects
        cpack_zhegvx_fallback.cpp
    )
  endif()

  if(USE_CUDA)
    set(cuda_objects ./kernels/hegvd_op.cpp ./kernels/cuda/diag_cusolver.cu
                     diago_cusolver.cpp diago_cusolver.h)
    if(ENABLE_CUSOLVERMP)
      list(APPEND cuda_objects ./kernels/cuda/diag_cusolvermp.cu
           diago_cusolvermp.cpp diago_cusolvermp.h)
    endif()
    add_library(diag_cusolver OBJECT ${cuda_objects})
    if(ENABLE_COVERAGE)
      add_coverage(diag_cusolver)
    endif()
  endif()

  if(ENABLE_PEXSI)
    list(APPEND objects diago_pexsi.cpp)
    add_subdirectory(module_pexsi)
  endif()

endif()

add_library(hsolver OBJECT ${objects})

option(ABACUS_DM_DIRECT_PILOT "Private cuSOLVER to Gamma FP64 DM trial" OFF)
if(ABACUS_DM_DIRECT_PILOT)
  if(NOT USE_CUDA OR NOT ENABLE_MPI OR NOT ENABLE_LCAO)
    message(FATAL_ERROR "ABACUS_DM_DIRECT_PILOT requires CUDA, MPI and LCAO")
  endif()
  target_compile_definitions(hsolver PRIVATE ABACUS_DM_DIRECT_PILOT)
  target_compile_definitions(diag_cusolver PRIVATE ABACUS_DM_DIRECT_PILOT)
endif()

if(ENABLE_COVERAGE)
  add_coverage(hsolver)
endif()

if(ENABLE_ELPA)
  add_subdirectory(module_genelpa)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
  if(ENABLE_MPI)
    add_subdirectory(kernels/test)
  endif()
  message(STATUS "Building tests")
endif()

add_subdirectory(module_linalg)
