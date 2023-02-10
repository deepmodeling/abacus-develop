find_program(BASH bash)
find_package(Python3 REQUIRED)
if(ENABLE_ASAN)
    add_test(
        NAME integrated_test_with_asan
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )

else()
    add_test(
        NAME integrated_test
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )
    add_executable(sum_ENV_H2 ${ABACUS_TEST_DIR}/integrate/tools/sum_ENV_H2.cpp)
    install(TARGETS sum_ENV_H2 DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)
    add_executable(sum_ENV_H2_cube ${ABACUS_TEST_DIR}/integrate/tools/sum_ENV_H2_cube.cpp)
    install(TARGETS sum_ENV_H2_cube DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)
    add_executable(sum_BAND_CHG_H2 ${ABACUS_TEST_DIR}/integrate/tools/sum_BAND_CHG_H2.cpp)
    install(TARGETS sum_BAND_CHG_H2 DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)
    add_executable(sum_BAND_CHG_H2_cube ${ABACUS_TEST_DIR}/integrate/tools/sum_BAND_CHG_H2_cube.cpp)
    install(TARGETS sum_BAND_CHG_H2_cube DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)

    if(ENABLE_COVERAGE)
      add_coverage(sum_ENV_H2)
      add_coverage(sum_ENV_H2_cube)
      add_coverage(sum_BAND_CHG_H2)
      add_coverage(sum_BAND_CHG_H2_cube)
    endif()
endif()
