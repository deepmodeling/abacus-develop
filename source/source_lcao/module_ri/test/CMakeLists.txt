abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES dm_mixing_test.cpp ../mix_dmk_2d.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES ri_cv_io_test.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_rpa_test
  LIBS ${math_libs}
  SOURCES sternheimer_rpa_test.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_hamiltonian_test
  LIBS ${math_libs}
  SOURCES sternheimer_fd_hamiltonian_test.cpp ../sternheimer_fd_hamiltonian.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_fd_solver_test
  LIBS ${math_libs}
  SOURCES sternheimer_fd_solver_test.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_dft_zero_order_test
  LIBS ${math_libs}
  SOURCES sternheimer_dft_zero_order_test.cpp ../sternheimer_dft_zero_order.cpp ../sternheimer_fd_solver.cpp ../sternheimer_fd_hamiltonian.cpp ../sternheimer_rpa.cpp
)
AddTest(
  TARGET MODULE_RI_sternheimer_abacus_fd_adapter_test
  LIBS base ${math_libs}
  SOURCES sternheimer_abacus_fd_adapter_test.cpp ../sternheimer_abacus_fd_adapter.cpp ../sternheimer_fd_hamiltonian.cpp
)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
