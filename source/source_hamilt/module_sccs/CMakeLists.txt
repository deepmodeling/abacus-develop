add_library(
  sccs OBJECT
  sccs_diag.cpp
  sccs_cavity.cpp
  sccs_cav_deriv.cpp
  sccs_parameters.cpp
  sccs_functional.cpp
  sccs_ion_charge.cpp
  sccs_lowpass.cpp
  sccs_ion_force.cpp
  sccs_coul_0d.cpp
  sccs_coul_2d.cpp
  sccs_moments.cpp
  sccs_pw_coulomb.cpp
  sccs_response.cpp)

if(ENABLE_COVERAGE)
  add_coverage(sccs)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
