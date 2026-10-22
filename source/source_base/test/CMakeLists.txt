abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(data)

addtest(
  TARGET MODULE_BASE_additional_coverage
  LIBS parameter base device container
  SOURCES additional_coverage_test.cpp)
addtest(
  TARGET MODULE_BASE_blas_connector
  LIBS parameter base device
  SOURCES blas_connector_test.cpp blas_connector_additional_test.cpp)
addtest(
  TARGET MODULE_BASE_atom_in
  LIBS parameter
  SOURCES atom_in_test.cpp)
addtest(
  TARGET MODULE_BASE_timer
  LIBS parameter
  SOURCES timer_test.cpp ../timer.cpp ../global_variable.cpp)
addtest(
  TARGET MODULE_BASE_tool_quit
  LIBS parameter
  SOURCES tool_quit_test.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../global_file.cpp
          ../global_function.cpp
          ../memory_recorder.cpp
          ../timer.cpp)
addtest(
  TARGET MODULE_BASE_tool_check
  LIBS parameter
  SOURCES tool_check_test.cpp
          ../tool_check.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../global_file.cpp
          ../global_function.cpp
          ../memory_recorder.cpp
          ../timer.cpp)
addtest(TARGET MODULE_BASE_memory
        SOURCES memory_test.cpp ../memory_recorder.cpp ../parallel_reduce.cpp)
addtest(
  TARGET MODULE_BASE_global_function
  LIBS parameter
  SOURCES global_function_test.cpp
          ../global_function.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../global_file.cpp
          ../memory_recorder.cpp
          ../timer.cpp)
addtest(
  TARGET MODULE_BASE_vector3
  LIBS parameter
  SOURCES vector3_test.cpp)
addtest(
  TARGET MODULE_BASE_matrix3
  LIBS parameter base device
  SOURCES matrix3_test.cpp)
addtest(
  TARGET MODULE_BASE_intarray
  LIBS parameter
  SOURCES intarray_test.cpp ../intarray.cpp)
addtest(
  TARGET MODULE_BASE_realarray
  LIBS parameter
  SOURCES realarray_test.cpp ../realarray.cpp)
addtest(
  TARGET MODULE_BASE_matrix
  LIBS parameter base device
  SOURCES matrix_test.cpp)
addtest(
  TARGET MODULE_BASE_complexarray
  LIBS parameter
  SOURCES complexarray_test.cpp
          ../complexarray.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../global_file.cpp
          ../global_function.cpp
          ../memory_recorder.cpp
          ../timer.cpp)
addtest(
  TARGET MODULE_BASE_complexmatrix
  LIBS parameter base device
  SOURCES complexmatrix_test.cpp)
addtest(
  TARGET MODULE_BASE_integral
  LIBS parameter
  SOURCES math_integral_test.cpp ../math_integral.cpp)
addtest(
  TARGET MODULE_BASE_sph_bessel_recursive
  LIBS parameter
  SOURCES sph_bessel_rec_test.cpp ../sph_bessel_rec_d1.cpp
          ../sph_bessel_rec_d2.cpp ../memory_recorder.cpp
          ../global_variable.cpp)
addtest(
  TARGET MODULE_BASE_ylmreal
  LIBS parameter base device
  SOURCES math_ylmreal_test.cpp ../libm/branred.cpp ../libm/sincos.cpp)
addtest(
  TARGET MODULE_BASE_math_sphbes
  LIBS parameter
  SOURCES math_sphbes_test.cpp ../math_sphbes.cpp ../timer.cpp)
addtest(
  TARGET MODULE_BASE_mathzone
  LIBS parameter base device
  SOURCES mathzone_test.cpp)
addtest(
  TARGET MODULE_BASE_mathzone_add1
  LIBS parameter base device
  SOURCES mathzone_add1_test.cpp)
addtest(
  TARGET MODULE_BASE_math_polyint
  LIBS parameter
  SOURCES math_polyint_test.cpp ../math_polyint.cpp ../realarray.cpp
          ../timer.cpp)
addtest(
  TARGET MODULE_BASE_gram_schmidt_orth
  LIBS parameter base device
  SOURCES gram_schmidt_orth_test.cpp)
addtest(
  TARGET MODULE_BASE_math_bspline
  LIBS parameter
  SOURCES math_bspline_test.cpp ../math_bspline.cpp)
addtest(
  TARGET MODULE_BASE_inverse_matrix
  LIBS parameter base device
  SOURCES inverse_matrix_test.cpp)
addtest(
  TARGET MODULE_BASE_mymath
  LIBS parameter
  SOURCES mymath_test.cpp ../mymath.cpp ../timer.cpp)
addtest(
  TARGET MODULE_BASE_container
  LIBS parameter
  SOURCES container_op_test.cpp ../container_operator.h)

addtest(
  TARGET MODULE_BASE_math_chebyshev
  LIBS parameter base device container
  SOURCES math_chebyshev_test.cpp)

addtest(
  TARGET MODULE_BASE_lapack_connector
  LIBS parameter base device
  SOURCES lapack_connector_test.cpp)

addtest(
  TARGET MODULE_BASE_opt_cg
  LIBS parameter base device
  SOURCES opt_cg_test.cpp opt_test_tools.cpp mpi_test_main.cpp)

addtest(TARGET MODULE_BASE_opt_dcsrch SOURCES test_opt_dcsrch.cpp
                                              ../opt_dcsrch.cpp)

addtest(
  TARGET MODULE_BASE_opt_tn
  LIBS parameter base device
  SOURCES opt_tn_test.cpp opt_test_tools.cpp mpi_test_main.cpp)

addtest(
  TARGET MODULE_BASE_ylm
  LIBS parameter
  SOURCES ylm_test.cpp
          ../ylm.cpp
          ../timer.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../global_file.cpp
          ../global_function.cpp
          ../memory_recorder.cpp)

addtest(
  TARGET MODULE_BASE_global_file
  LIBS parameter
  SOURCES global_file_test.cpp
          ../global_file.cpp
          ../global_function.cpp
          ../tool_quit.cpp
          ../global_variable.cpp
          ../memory_recorder.cpp
          ../timer.cpp)

addtest(
  TARGET MODULE_BASE_tool_title
  LIBS parameter
  SOURCES tool_title_test.cpp
          ../tool_title.cpp
          ../global_variable.cpp
          ../global_function.cpp
          ../timer.cpp
          ../tool_quit.cpp
          ../global_file.cpp
          ../memory_recorder.cpp)

addtest(
  TARGET MODULE_BASE_element_basis_index
  LIBS parameter
  SOURCES elem_basis_idx_test.cpp ../element_basis_index.cpp)

addtest(
  TARGET MODULE_BASE_tool_threading
  LIBS parameter
  SOURCES tool_threading_test.cpp ../tool_threading.h)

addtest(
  TARGET MODULE_BASE_sph_bessel_tf
  SOURCES sph_bessel_tf_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_cubic_spline
  SOURCES cubic_spline_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_math_erf_complex
  SOURCES math_erf_complex_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_math_lib_info
  SOURCES math_lib_info_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_clebsch_gordan_coeff_test
  SOURCES cg_coeff_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_assoc_laguerre_test
  SOURCES assoc_laguerre_test.cpp
  LIBS parameter base device)

addtest(
  TARGET MODULE_BASE_ndarray_test
  LIBS parameter
  SOURCES ndarray_test.cpp)

addtest(
  TARGET MODULE_BASE_formatter_test
  LIBS parameter
  SOURCES formatter_test.cpp)

addtest(TARGET MODULE_BASE_lebedev_laikov
        SOURCES test_lebedev_laikov.cpp ../ylm.cpp ../math_lebedev_laikov.cpp)

abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_binstream
  LIBS parameter base device
  SOURCES binstream_test.cpp)

addtest(
  TARGET MODULE_BASE_sparse_matrix
  LIBS parameter base device
  SOURCES sparse_matrix_test.cpp)

addtest(
  TARGET MODULE_BASE_file_reader
  LIBS parameter base device
  SOURCES file_reader_test.cpp)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES csr_reader_test.cpp)

if(ENABLE_GOOGLEBENCH)
  addtest(
    TARGET MODULE_BASE_perf_sphbes
    LIBS parameter
    SOURCES perf_sphbes_test.cpp ../math_sphbes.cpp ../timer.cpp)
endif()
