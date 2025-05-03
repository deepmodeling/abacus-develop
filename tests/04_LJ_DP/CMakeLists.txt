enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 04_LJ_DP_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/04_LJ_DP
    )
else()
    add_test(
        NAME 04_LJ_DP
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/04_LJ_DP
    )
endif()
