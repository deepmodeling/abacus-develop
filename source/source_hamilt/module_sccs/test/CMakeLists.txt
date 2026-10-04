AddTest(TARGET MODULE_HAMILT_sccs_cavity
  SOURCES test_sccs_cavity.cpp ../sccs_cavity.cpp)
AddTest(TARGET MODULE_HAMILT_sccs_parameters
  SOURCES test_sccs_parameters.cpp ../sccs_parameters.cpp ../sccs_cavity.cpp)
