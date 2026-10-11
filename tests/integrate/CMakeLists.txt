find_program(BASH bash)
find_package(Python3 REQUIRED)

add_test(
    NAME compare_file_test
    COMMAND ${Python3_EXECUTABLE} -m unittest tools/test_compare_file.py
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
)

add_test(
    NAME compare_wannier90_test
    COMMAND ${Python3_EXECUTABLE} -m unittest tools/test_compare_wannier90.py
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
)

add_test(
    NAME catch_properties_test
    COMMAND ${BASH} tools/test_catch_properties.sh
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/integrate
)
   
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
