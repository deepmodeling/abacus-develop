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

# Autotest checks numerical results; this existing driver also asserts the
# structure-file inventory for independent output and restart cadences.
find_program(ABACUS_CADENCE_PYTHON python3)
if(ABACUS_CADENCE_PYTHON AND NOT ENABLE_ASAN)
  add_test(
    NAME 04_FF_md_output_cadence
    COMMAND
      ${ABACUS_CADENCE_PYTHON}
      ${ABACUS_TEST_DIR}/integrate/validation_tools/check_md_output_cadence.py
      ${ABACUS_BIN_PATH} ${ABACUS_TEST_DIR}/04_FF/21_LJ_cadence)
  set_tests_properties(04_FF_md_output_cadence
                       PROPERTIES TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1")
endif()
