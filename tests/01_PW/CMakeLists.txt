enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 01_PW_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW
    )
else()
    add_test(
        NAME 01_PW
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW
    )
endif()
