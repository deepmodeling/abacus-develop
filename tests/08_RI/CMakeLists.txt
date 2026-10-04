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

    # These compact cases validate only producer artifacts consumed by LibRPA.
    # They leave the established EXX and BSE regression corpus unchanged.
    add_test(
        NAME 08_RI_librpa_producer
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 1 -j auto -f CASES_LIBRPA_PRODUCER.txt
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI
    )
endif()
