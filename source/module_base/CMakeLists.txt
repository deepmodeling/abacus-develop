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
    blas_connector.cpp
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
    math_integral.cpp
    math_lebedev_laikov.cpp
    math_polyint.cpp
    math_sphbes.cpp
    math_ylmreal.cpp
    math_bspline.cpp
    math_chebyshev.cpp
    mathzone_add1.cpp
    matrix.cpp
    matrix3.cpp
    memory.cpp
    mymath.cpp
    opt_CG.cpp
    opt_DCsrch.cpp
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
    module_mixing/mixing_data.cpp
    module_mixing/mixing.cpp
    module_mixing/plain_mixing.cpp
    module_mixing/pulay_mixing.cpp
    module_mixing/broyden_mixing.cpp
    ${LIBM_SRC}
)

target_link_libraries(base PUBLIC container)

add_subdirectory(module_container)

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
    add_subdirectory(grid/test)
    if (USE_ABACUS_LIBM)
      add_subdirectory(libm/test)
    endif()
  endif()
endif()
