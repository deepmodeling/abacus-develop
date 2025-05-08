enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 13_NAO_multik_GPU_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -f CASES_GPU.txt -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/13_NAO_multik_GPU
    )
else()
    add_test(
        NAME 13_NAO_multik_GPU
        COMMAND ${BASH} ../integrate/Autotest.sh -f CASES_GPU.txt -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/13_NAO_multik_GPU
    )
endif()
