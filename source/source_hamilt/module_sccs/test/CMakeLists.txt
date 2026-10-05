AddTest(TARGET MODULE_HAMILT_sccs_cavity
  SOURCES test_sccs_cavity.cpp ../sccs_cavity.cpp)
AddTest(TARGET MODULE_HAMILT_sccs_parameters
  SOURCES test_sccs_parameters.cpp ../sccs_parameters.cpp ../sccs_cavity.cpp)

AddTest(TARGET MODULE_HAMILT_sccs_poisson
  LIBS parameter planewave device base
  SOURCES test_sccs_pw_coulomb.cpp test_sccs_ionic_charge.cpp test_sccs_ionic_force.cpp test_sccs_response.cpp test_sccs_functional.cpp
    ../../../source_estate/test/test_pot_sccs.cpp
    ../../../source_estate/module_pot/pot_sccs.cpp
    ../../../source_cell/cell_tools.cpp ../../../source_cell/atom_spec.cpp
    ../../../source_cell/atom_pseudo.cpp ../../../source_cell/pseudo.cpp
    ../sccs_functional.cpp
    ../sccs_parameters.cpp
    ../sccs_ionic_charge.cpp ../sccs_ionic_force.cpp ../sccs_pw_coulomb.cpp ../sccs_response.cpp ../sccs_cavity.cpp
    ../../module_xc/xc_grad_utils.cpp)
set_tests_properties(MODULE_HAMILT_sccs_poisson PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60)

if(ENABLE_MPI)
  add_test(NAME MODULE_HAMILT_sccs_poisson_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:MODULE_HAMILT_sccs_poisson> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_HAMILT_sccs_poisson_mpi PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 60 PROCESSORS 2)
endif()
