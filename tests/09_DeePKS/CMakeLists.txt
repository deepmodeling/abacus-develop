enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 09_DeePKS_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/09_DeePKS
    )
    add_test(
        NAME 09_DeePKS_test1_with_asan
        COMMAND ${BASH} Autotest1.sh -a ${CMAKE_CURRENT_BINARY_DIR}/../../source/module_hamilt_lcao/module_deepks/test/test_deepks
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/09_DeePKS
    )

else()
    add_test(
        NAME 09_DeePKS
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/09_DeePKS
    )
    # TODO: I will rewrite the unit tests and remove 604 to module_deepks/test/
    add_test(
       NAME 09_DeePKS_test1
       COMMAND ${BASH} Autotest1.sh -a ${CMAKE_CURRENT_BINARY_DIR}/../../source/module_hamilt_lcao/module_deepks/test/test_deepks
       WORKING_DIRECTORY ${ABACUS_TEST_DIR}/09_DeePKS
    )
endif()
