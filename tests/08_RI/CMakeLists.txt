enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 08_RI_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI
    )
else()
    add_test(
        NAME 08_RI
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI
    )
endif()
