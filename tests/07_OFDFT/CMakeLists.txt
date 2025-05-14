enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 07_OFDFT_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/07_OFDFT
    )
else()
    add_test(
        NAME 07_OFDFT
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/07_OFDFT
    )
endif()
