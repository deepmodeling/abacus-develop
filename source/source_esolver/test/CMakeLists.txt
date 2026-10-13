abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__LCAO)

AddTest(
  TARGET MODULE_ESOLVER_electronic_occupations_test
  LIBS base device MPI::MPI_CXX
  KEEP_FEATURE_DEFINITIONS __MPI
  SOURCES test_electronic_occupations.cpp ../electronic_occupations.cpp
)
set_tests_properties(MODULE_ESOLVER_electronic_occupations_test
  PROPERTIES TIMEOUT 30 ENVIRONMENT "OMP_NUM_THREADS=1")
add_test(
  NAME MODULE_ESOLVER_electronic_occupations_mpi4_test
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_ESOLVER_electronic_occupations_test> ${MPIEXEC_POSTFLAGS}
)
set_tests_properties(MODULE_ESOLVER_electronic_occupations_mpi4_test
  PROPERTIES TIMEOUT 30 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=1")

set(_esolver_dp_test_libs
  parameter
  base
  device)

if(DEFINED DeePMD_DIR)
  if(DeePMDC_FOUND)
    list(APPEND _esolver_dp_test_libs DeePMD::deepmd_c)
  else()
    list(APPEND _esolver_dp_test_libs DeePMD::deepmd_cc)
  endif()
endif()

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_ESOLVER_esolver_dp_test
  LIBS ${_esolver_dp_test_libs}
  SOURCES
    esolver_dp_test.cpp
    ../esolver_dp.cpp
    ../../source_cell/basecell.cpp
    ../../source_cell/cif_io.cpp
    ../../source_cell/module_neighlist/neighbor_list.cpp
    ../../source_io/module_output/output_log.cpp
)
