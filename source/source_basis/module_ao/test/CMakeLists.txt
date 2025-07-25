remove_definitions(-D__EXX)

list(APPEND depend_files
 ../../../source_base/math_integral.cpp
  ../../../source_base/math_sphbes.cpp
  ../../../source_base/math_polyint.cpp
  ../../../source_base/math_ylmreal.cpp
  ../../../source_base/ylm.cpp
  ../../../source_base/memory.cpp
  ../../../source_base/module_external/blas_connector_base.cpp
  ../../../source_base/module_external/blas_connector_vector.cpp
  ../../../source_base/module_external/blas_connector_matrix.cpp
  ../../../source_base/complexarray.cpp
  ../../../source_base/complexmatrix.cpp
  ../../../source_base/matrix.cpp
  ../../../source_base/realarray.cpp
  ../../../source_base/intarray.cpp
  ../../../source_base/sph_bessel_recursive-d1.cpp
  ../../../source_base/sph_bessel_recursive-d2.cpp
  ../../../source_base/tool_title.cpp
  ../../../source_base/tool_quit.cpp
  ../../../source_base/tool_check.cpp
  ../../../source_base/timer.cpp
  ../../../source_base/mathzone_add1.cpp
  ../../../source_base/global_variable.cpp
  ../../../source_base/global_function.cpp
  ../../../source_base/global_file.cpp
  ../../../source_base/libm/branred.cpp
  ../../../source_base/libm/sincos.cpp
  ../../../source_base/spherical_bessel_transformer.cpp
  ../../../source_base/cubic_spline.cpp
  ../../../source_base/parallel_2d.cpp
  ../ORB_read.cpp
  ../ORB_atomic.cpp
  ../ORB_atomic_lm.cpp
  ../ORB_nonlocal.cpp
  ../ORB_nonlocal_lm.cpp
  ../ORB_gaunt_table.cpp
  ../parallel_orbitals.cpp
  ../../../source_lcao/center2_orb-orb11.cpp
  ../../../source_lcao/center2_orb.cpp
  )

install(DIRECTORY GaAs DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../../tests)
install(DIRECTORY GaAs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


AddTest(
    TARGET ORB_nonlocal_test
    SOURCES ORB_nonlocal_test.cpp
        ../ORB_nonlocal.cpp
        ../ORB_nonlocal_lm.cpp
    LIBS parameter ${math_libs} device base
)

AddTest(
    TARGET ORB_nonlocal_lm_test
    SOURCES ORB_nonlocal_lm_test.cpp
        ../ORB_nonlocal_lm.cpp
    LIBS parameter ${math_libs} device base
)

AddTest(
    TARGET ORB_atomic_lm_test
    SOURCES ORB_atomic_lm_test.cpp
        ../ORB_atomic_lm.cpp
    LIBS parameter ${math_libs} device base
)

AddTest(
    TARGET ORB_read_test
    SOURCES ORB_read_test.cpp
        ../ORB_read.cpp
        ../ORB_atomic.cpp
        ../ORB_atomic_lm.cpp
    LIBS parameter ${math_libs} device base 
)


AddTest(
    TARGET parallel_orbitals_test 
    SOURCES parallel_orbitals_test.cpp ../parallel_orbitals.cpp
    LIBS parameter ${math_libs} device base
)

install(FILES parallel_orbitals_test.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
find_program(BASH bash)
add_test(NAME parallel_orbitals_test_para
      COMMAND ${BASH} parallel_orbitals_test.sh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

install(DIRECTORY lcao_H2O DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(DIRECTORY lcao_H2O DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../../tests)

