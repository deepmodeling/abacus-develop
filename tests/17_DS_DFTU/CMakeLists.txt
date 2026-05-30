enable_testing()

# TODO: Enable all tests after DeltaSpin + DFT+U code is merged into develop.
# These tests require the DeltaSpin module (source/module_hamilt_general/hamilt_deltaspin.h)
# and associated DFT+U with spin constraint functionality.
# Uncomment the add_test blocks below once the feature branch is merged.

find_program(BASH bash)
if(ENABLE_ASAN)
    add_test(
        NAME 17_DS_DFTU_test_with_asan
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_DS_DFTU
    )
else()
    add_test(
        NAME 17_DS_DFTU
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_DS_DFTU
    )
endif()
