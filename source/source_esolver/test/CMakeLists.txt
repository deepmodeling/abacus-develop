abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__LCAO)

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
    ../../source_io/module_output/cif_io.cpp
    ../../source_io/module_output/output_log.cpp
)
