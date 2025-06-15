enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 10_others_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/10_others
    )
else()
    add_test(
        NAME 10_others
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/10_others
    )
endif()
