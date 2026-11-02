add_library(
  module_rt_pw OBJECT
  cn_solver_pw.cpp
  td_orbital_ops.cpp
  td_prop_pw.cpp
  td_hamilt_pw.cpp
  prop_cn_pw.cpp
  prop_etrs_pw.cpp
  td_history_pw.cpp
  td_conserv_pw.cpp
  td_endpoint_pw.cpp
  td_pw.cpp
  td_vel_pw.cpp
  td_spin_pw.cpp
  op_pw_ekin_td.cpp
  op_pw_nl_td.cpp)

if(ENABLE_COVERAGE)
  add_coverage(module_rt_pw)
endif()

if(BUILD_TESTING)
  add_subdirectory(unittests)
endif()
