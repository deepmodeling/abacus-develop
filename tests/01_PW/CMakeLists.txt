enable_testing()

find_program(BASH bash)
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

# rVV10 with a built-in functional
if(ENABLE_LIBVDWXC)
  add_test(
    NAME 01_PW_libvdwxc
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
            -f CASES_LIBVDWXC.txt
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW)
endif()

# rVV10 with the Libxc functional rPW86-PBE
if(ENABLE_LIBVDWXC AND ENABLE_LIBXC)
  add_test(
    NAME 01_PW_libvdwxc_libxc
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
            -f CASES_LIBVDWXC_LIBXC.txt
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/01_PW)
endif()
