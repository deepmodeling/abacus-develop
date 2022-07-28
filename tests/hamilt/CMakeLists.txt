add_test(
    NAME hamilt_test_kinetic
    COMMAND ${BASH} Autotest.sh -a ${CMAKE_CURRENT_BINARY_DIR}/../../source/module_hamilt/test/hamilt_test_kinetic -r ${ABACUS_TEST_DIR}/hamilt/kin
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/hamilt
)