if (USE_ABACUS_LIBM)
list (APPEND LIBM_SRC
  libm/branred.cpp
  libm/cexp.cpp
  libm/exp.cpp
  libm/sincos.cpp
)
endif()
add_library(
    base
    OBJECT
    assoc_laguerre.cpp
    module_external/blas_connector_base.cpp
    module_external/blas_connector_vector.cpp
    module_external/blas_connector_matrix.cpp
    clebsch_gordan_coeff.cpp
    complexarray.cpp
    complexmatrix.cpp
    element_basis_index.cpp
    export.cpp
    inverse_matrix.cpp
    gather_math_lib_info.cpp
    global_file.cpp
    global_function.cpp
    global_function_ddotreal.cpp
    global_variable.cpp
    intarray.cpp
    mathzone_add1.cpp
    matrix.cpp
    matrix3.cpp
    memory.cpp
    opt_CG.cpp
    opt_DCsrch.cpp
    para_gemm.cpp
    realarray.cpp
    sph_bessel_recursive-d1.cpp
    sph_bessel_recursive-d2.cpp
    timer.cpp
    tool_check.cpp
    tool_quit.cpp
    tool_title.cpp
    ylm.cpp
    abfs-vector3_order.cpp
    parallel_common.cpp
    parallel_global.cpp
    parallel_comm.cpp
    parallel_reduce.cpp
    parallel_device.cpp
    spherical_bessel_transformer.cpp
    cubic_spline.cpp
    parallel_2d.cpp
    projgen.cpp
    module_mixing/mixing_data.cpp
    module_mixing/mixing.cpp
    module_mixing/plain_mixing.cpp
    module_mixing/pulay_mixing.cpp
    module_mixing/broyden_mixing.cpp
    module_math/math_bspline.cpp
    module_math/math_chebyshev.cpp
    module_math/math_erf_complex.cpp
    module_math/math_integral.cpp
    module_math/math_lebedev_laikov.cpp
    module_math/math_polyint.cpp
    module_math/math_sphbes.cpp
    module_math/math_ylmreal.cpp
    module_math/mymath.cpp
    ${LIBM_SRC}
)

target_link_libraries(base PUBLIC container device)

target_include_directories(base
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/module_math
)
if (USE_DSP)
  target_link_libraries(base PUBLIC ${MTBLAS_FFT_DIR}/libmtblas/lib/libmtblas.a)
  target_link_libraries(base PUBLIC ${MTBLAS_FFT_DIR}/libmtblas/lib/libmtblasdev.a)
endif()
add_subdirectory(module_container)
add_subdirectory(module_math)

if(ENABLE_COVERAGE)
  add_coverage(base)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI) #Zhang Xiaoyang: A temporary solution. Please classify the serial and parrallel tests later. 2024-7-8
    add_subdirectory(test)
    add_subdirectory(test_parallel)
    add_subdirectory(kernels/test)
    add_subdirectory(module_mixing/test)
    add_subdirectory(module_device/test)
    add_subdirectory(module_grid/test)
    if (USE_ABACUS_LIBM)
      add_subdirectory(libm/test)
    endif()
  endif()
endif()
