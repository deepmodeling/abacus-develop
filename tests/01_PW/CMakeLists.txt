enable_testing()

find_program(BASH bash)
if(ENABLE_MPI)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(
    NAME 01_PW_dftu_first_step_restart
    COMMAND ${Python3_EXECUTABLE}
            ${ABACUS_TEST_DIR}/integrate/validation_tools/test_dftu_first_step_restart.py
            ${ABACUS_BIN_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/relax_dftu_spin2)
endif()
if(ENABLE_ASAN)
  add_test(
    NAME 01_PW_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW)
else()
  add_test(
    NAME 01_PW
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW)
endif()
