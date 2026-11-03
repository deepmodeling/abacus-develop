enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
  add_test(
    NAME 04_FF_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/04_FF)
else()
  add_test(
    NAME 04_FF
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/04_FF)
endif()

find_program(SNAP_TEST_PYTHON python3)
if(ENABLE_MPI AND SNAP_TEST_PYTHON AND MPIEXEC_EXECUTABLE)
  add_test(
    NAME 04_FF_SNAP_rank_error
    COMMAND ${SNAP_TEST_PYTHON} ${ABACUS_TEST_DIR}/integrate/validation_tools/check_snap_rank_error.py
            --abacus ${ABACUS_BIN_PATH}
            --mpiexec ${MPIEXEC_EXECUTABLE}
            --numproc-flag=${MPIEXEC_NUMPROC_FLAG}
            --models ${ABACUS_TEST_DIR}/SIMULATION_MODELS)
  set_tests_properties(04_FF_SNAP_rank_error PROPERTIES TIMEOUT 120)
endif()
