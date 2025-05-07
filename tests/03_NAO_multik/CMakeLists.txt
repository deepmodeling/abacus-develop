enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 03_NAO_multik_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik
    )
else()
    add_test(
        NAME 03_NAO_multik
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik
    )
endif()
