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
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(
        NAME 08_RI_librpa_wavefunction_checks
        COMMAND ${Python3_EXECUTABLE} -m unittest discover
                -s ${ABACUS_TEST_DIR}/integrate/tools -p test_librpa_wavefunctions.py
    )
    add_test(
        NAME 08_RI_librpa_producer_mpi
        COMMAND ${Python3_EXECUTABLE} ${ABACUS_TEST_DIR}/integrate/tools/run_librpa_mpi.py
                --abacus ${ABACUS_BIN_PATH} --cases ${ABACUS_TEST_DIR}/08_RI
                --mpirun ${MPIEXEC_EXECUTABLE} --mpi-np-flag ${MPIEXEC_NUMPROC_FLAG}
    )
    set_tests_properties(08_RI_librpa_producer_mpi PROPERTIES PROCESSORS 2 TIMEOUT 1200)
    add_test(
        NAME 08_RI_librpa_symmetry
        COMMAND ${Python3_EXECUTABLE} ${ABACUS_TEST_DIR}/integrate/tools/run_librpa_sym.py
                --abacus ${ABACUS_BIN_PATH} --cases ${ABACUS_TEST_DIR}/08_RI
                --mpirun ${MPIEXEC_EXECUTABLE} --mpi-np-flag ${MPIEXEC_NUMPROC_FLAG}
    )
    set_tests_properties(08_RI_librpa_symmetry PROPERTIES PROCESSORS 2 TIMEOUT 1800)
endif()
