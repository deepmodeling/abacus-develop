enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
	NAME 17_TDOFDFT_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
	WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_TDOFDFT
    )
else()
    add_test(
	NAME 17_TDOFDFT
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
	WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_TDOFDFT
    )
endif()
