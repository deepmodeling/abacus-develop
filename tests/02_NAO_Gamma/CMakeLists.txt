enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 02_NAO_Gamma_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/02_NAO_Gamma
    )
else()
    add_test(
        NAME 02_NAO_Gamma
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/02_NAO_Gamma
    )
endif()
