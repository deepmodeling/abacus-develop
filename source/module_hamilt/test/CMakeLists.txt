remove_definitions(-D__MPI)
remove_definitions(-D__EXX)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__DEEPKS)
remove_definitions(-D__LCAO)

add_executable(
  hamilt_test_kinetic
  kinetic_test.cpp ../ks_pw/ekinetic_pw.cpp ../../module_psi/psi.cpp
)

target_link_libraries(
  hamilt_test_kinetic
  ${math_libs} planewave_serial base_serial
)

install(
  TARGETS hamilt_test_kinetic
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../../tests
)