enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 08_EXX_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_EXX
    )
else()
    add_test(
        NAME 08_EXX
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_EXX
    )
    # The SOC-HSE case uses substantial per-thread memory; retain four MPI ranks
    # while avoiding an OOM kill on the CI runner.
    set_tests_properties(08_EXX PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
endif()
