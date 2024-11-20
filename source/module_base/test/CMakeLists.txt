remove_definitions(-D__MPI)
install(DIRECTORY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
AddTest(
  TARGET base_blas_connector
  LIBS parameter  ${math_libs}
  SOURCES blas_connector_test.cpp ../blas_connector.cpp
)
AddTest(
  TARGET base_atom_in
  LIBS parameter 
  SOURCES atom_in_test.cpp
)
AddTest(
  TARGET base_timer
  LIBS parameter 
  SOURCES timer_test.cpp ../timer.cpp  ../global_variable.cpp
)
AddTest(
  TARGET base_tool_quit
  LIBS parameter 
  SOURCES tool_quit_test.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp  
)
AddTest(
  TARGET base_tool_check
  LIBS parameter 
  SOURCES tool_check_test.cpp ../tool_check.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp 
)
AddTest(
  TARGET base_memory
  SOURCES memory_test.cpp ../memory.cpp ../parallel_reduce.cpp
)
ADDTest(
  TARGET base_global_function
  LIBS parameter ${math_libs} 
  SOURCES global_function_test.cpp ../blas_connector.cpp ../global_function.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_vector3
  LIBS parameter 
  SOURCES vector3_test.cpp
)
AddTest(
  TARGET base_matrix3
  LIBS parameter  ${math_libs}
  SOURCES matrix3_test.cpp ../matrix3.cpp ../matrix.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp ../blas_connector.cpp
)
AddTest(
  TARGET base_intarray
  LIBS parameter 
  SOURCES intarray_test.cpp ../intarray.cpp
)
AddTest(
  TARGET base_realarray
  LIBS parameter 
  SOURCES realarray_test.cpp ../realarray.cpp
)
AddTest(
 TARGET base_matrix
  LIBS parameter  ${math_libs}
  SOURCES matrix_test.cpp ../blas_connector.cpp ../matrix.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_complexarray
  LIBS parameter 
  SOURCES complexarray_test.cpp ../complexarray.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp
)
AddTest(
  TARGET base_complexmatrix
  LIBS parameter  ${math_libs}
  SOURCES complexmatrix_test.cpp ../blas_connector.cpp ../complexmatrix.cpp ../matrix.cpp
)
AddTest(
  TARGET base_integral
  LIBS parameter 
  SOURCES math_integral_test.cpp ../math_integral.cpp
)
AddTest(
  TARGET base_sph_bessel_recursive
  LIBS parameter 
  SOURCES sph_bessel_recursive_test.cpp ../sph_bessel_recursive-d1.cpp ../sph_bessel_recursive-d2.cpp
)
AddTest(
  TARGET base_ylmreal
  LIBS parameter  ${math_libs} device
  SOURCES math_ylmreal_test.cpp ../blas_connector.cpp ../math_ylmreal.cpp ../complexmatrix.cpp ../global_variable.cpp ../ylm.cpp ../realarray.cpp ../timer.cpp ../matrix.cpp ../vector3.h
          ../parallel_reduce.cpp ../parallel_global.cpp ../parallel_comm.cpp ../parallel_common.cpp
          ../memory.cpp ../libm/branred.cpp ../libm/sincos.cpp
)
AddTest(
  TARGET base_math_sphbes
  LIBS parameter 
  SOURCES math_sphbes_test.cpp ../math_sphbes.cpp ../timer.cpp
)
AddTest(
  TARGET base_mathzone
  LIBS parameter  ${math_libs}
  SOURCES mathzone_test.cpp ../matrix3.cpp ../matrix.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp ../blas_connector.cpp
)
AddTest(
  TARGET base_mathzone_add1
  LIBS parameter  ${math_libs}
  SOURCES mathzone_add1_test.cpp ../blas_connector.cpp ../mathzone_add1.cpp ../math_sphbes.cpp ../matrix3.cpp ../matrix.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp 
)
AddTest(
  TARGET base_math_polyint
  LIBS parameter 
  SOURCES math_polyint_test.cpp ../math_polyint.cpp ../realarray.cpp ../timer.cpp
)
AddTest(
  TARGET base_gram_schmidt_orth
  LIBS parameter  ${math_libs}
  SOURCES gram_schmidt_orth_test.cpp ../blas_connector.cpp ../gram_schmidt_orth.h ../gram_schmidt_orth-inl.h ../global_function.h ../math_integral.cpp
)
AddTest(
  TARGET base_math_bspline
  LIBS parameter 
  SOURCES math_bspline_test.cpp ../math_bspline.cpp
)
AddTest(
  TARGET base_inverse_matrix
  LIBS parameter  ${math_libs}
  SOURCES inverse_matrix_test.cpp ../blas_connector.cpp ../inverse_matrix.cpp ../complexmatrix.cpp ../matrix.cpp ../timer.cpp 
)
AddTest(
  TARGET base_mymath
  LIBS parameter 
  SOURCES mymath_test.cpp ../mymath.cpp ../timer.cpp
)
AddTest(
  TARGET base_container
  LIBS parameter 
  SOURCES container_operator_test.cpp ../container_operator.h
)

AddTest(
  TARGET base_math_chebyshev
  LIBS parameter  ${math_libs} device container
  SOURCES math_chebyshev_test.cpp  ../blas_connector.cpp ../math_chebyshev.cpp ../tool_quit.cpp ../global_variable.cpp ../timer.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../parallel_reduce.cpp
)

AddTest(
  TARGET base_lapack_connector
  LIBS parameter  ${math_libs}
  SOURCES lapack_connector_test.cpp ../blas_connector.cpp ../lapack_connector.h
)

AddTest(
  TARGET base_opt_CG
  LIBS parameter  ${math_libs}
  SOURCES opt_CG_test.cpp opt_test_tools.cpp ../blas_connector.cpp ../opt_CG.cpp ../opt_DCsrch.cpp ../global_variable.cpp ../parallel_reduce.cpp
)

AddTest(
  TARGET base_opt_TN
  LIBS parameter  ${math_libs}
  SOURCES opt_TN_test.cpp opt_test_tools.cpp ../blas_connector.cpp ../opt_CG.cpp ../opt_DCsrch.cpp ../global_variable.cpp ../parallel_reduce.cpp
)

AddTest(
  TARGET base_ylm
  LIBS parameter 
  SOURCES ylm_test.cpp  ../ylm.cpp  ../timer.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp
)

AddTest(
  TARGET base_abfs-vector3_order
  LIBS parameter 
  SOURCES abfs-vector3_order_test.cpp ../abfs-vector3_order.cpp
)

AddTest(
  TARGET base_global_file
  LIBS parameter 
  SOURCES global_file_test.cpp ../global_file.cpp ../global_function.cpp ../tool_quit.cpp ../global_variable.cpp ../memory.cpp ../timer.cpp
)

AddTest(
  TARGET base_tool_title
  LIBS parameter 
  SOURCES tool_title_test.cpp ../tool_title.cpp  ../global_variable.cpp ../global_function.cpp ../timer.cpp  ../tool_quit.cpp  ../global_file.cpp ../memory.cpp
)

AddTest(
  TARGET base_element_basis_index
  LIBS parameter 
  SOURCES element_basis_index_test.cpp ../element_basis_index.cpp
)

AddTest(
  TARGET base_tool_threading
  LIBS parameter 
  SOURCES tool_threading_test.cpp ../tool_threading.h
)

AddTest(
  TARGET spherical_bessel_transformer
  SOURCES spherical_bessel_transformer_test.cpp ../blas_connector.cpp ../spherical_bessel_transformer.cpp ../math_sphbes.cpp ../math_integral.cpp ../timer.cpp
  LIBS parameter  ${math_libs}
)

AddTest(
  TARGET cubic_spline
  SOURCES cubic_spline_test.cpp ../blas_connector.cpp ../cubic_spline.cpp
  LIBS parameter  ${math_libs}
)

AddTest(
  TARGET clebsch_gordan_coeff_test
  SOURCES clebsch_gordan_coeff_test.cpp ../blas_connector.cpp ../clebsch_gordan_coeff.cpp ../intarray.cpp ../realarray.cpp ../complexmatrix.cpp ../matrix.cpp ../timer.cpp
  ../math_ylmreal.cpp  ../global_variable.cpp ../ylm.cpp  ../timer.cpp  ../vector3.h ../parallel_reduce.cpp ../parallel_global.cpp ../parallel_comm.cpp ../parallel_common.cpp
  ../memory.cpp ../libm/branred.cpp ../libm/sincos.cpp ../inverse_matrix.cpp ../lapack_connector.h
  LIBS parameter  ${math_libs} device
)

AddTest(
  TARGET assoc_laguerre_test
  SOURCES assoc_laguerre_test.cpp ../blas_connector.cpp ../assoc_laguerre.cpp ../tool_quit.cpp ../global_variable.cpp ../global_file.cpp ../global_function.cpp ../memory.cpp ../timer.cpp
  LIBS parameter  ${math_libs}
)

AddTest(
  TARGET ndarray_test
  LIBS parameter 
  SOURCES ndarray_test.cpp
)

AddTest(
  TARGET formatter_test
  LIBS parameter 
  SOURCES formatter_test.cpp
)

AddTest(
  TARGET lebedev_laikov
  SOURCES test_lebedev_laikov.cpp ../ylm.cpp ../math_lebedev_laikov.cpp
)

if(ENABLE_GOOGLEBENCH)
  AddTest(
    TARGET perf_sphbes
  LIBS parameter 
    SOURCES perf_sphbes_test.cpp ../math_sphbes.cpp ../timer.cpp
  )
endif()
