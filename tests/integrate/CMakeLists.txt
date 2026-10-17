find_program(BASH bash)
# Validation-script regressions do not impose a Python configure dependency
# on builds that only use a dedicated output validator.
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(
        NAME compare_file_test
        COMMAND ${Python3_EXECUTABLE} -m unittest validation_tools/test_compare_file.py
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )

    add_test(
        NAME catch_properties_test
        COMMAND ${BASH} validation_tools/test_catch_properties.sh
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )
else()
    message(STATUS "Python3 interpreter not found; skipping validation-script regression tests")
endif()
   
if(ENABLE_ASAN)
    add_test(
        NAME integrated_test_with_asan
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )

else()
    add_test(
        NAME integrated_test
        COMMAND ${BASH} Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
    )
endif()
