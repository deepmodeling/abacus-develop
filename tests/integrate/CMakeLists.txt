find_program(BASH bash)
find_package(Python3 REQUIRED)
   
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
if(ENABLE_MPI AND ENABLE_LIBXC)
    add_test(
        NAME PW_EXX_UNEVEN_KPAR
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_pw_exx_kpar.py
                --binary ${ABACUS_BIN_PATH}
                --source ${ABACUS_TEST_DIR}
                --launcher ${MPIEXEC_EXECUTABLE}
    )
    set_tests_properties(PW_EXX_UNEVEN_KPAR PROPERTIES TIMEOUT 2400 PROCESSORS 6)
endif()
