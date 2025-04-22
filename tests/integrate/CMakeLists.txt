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
    add_executable(sum_cube.exe ${ABACUS_TEST_DIR}/integrate/tools/sum_cube.cpp)
    install(TARGETS sum_cube.exe DESTINATION ${ABACUS_TEST_DIR}/integrate/tools/)

    if(ENABLE_COVERAGE)
      add_coverage(sum_cube.exe)
    endif()
endif()
