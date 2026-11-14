enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
  add_test(
    NAME 03_NAO_multik_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
else()
  add_test(
    NAME 03_NAO_multik
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
endif()

if(ENABLE_CNPY)
  add_test(
    NAME 03_NAO_multik_npz
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
            -f CASES_CNPY.txt
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
endif()

# rVV10 with the Libxc functional rPW86-PBE
if(ENABLE_LIBVDWXC AND ENABLE_LIBXC)
  add_test(
    NAME 03_NAO_multik_libvdwxc_libxc
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
            -f CASES_LIBVDWXC_LIBXC.txt
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
endif()
