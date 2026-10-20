find_program(BASH bash)
   
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

# Cheap collector regressions; no ABACUS calculation or finite differences.
find_program(LR_FORCE_TEST_PYTHON python3)
if(LR_FORCE_TEST_PYTHON)
    add_test(
        NAME lr_force_collector
        COMMAND ${LR_FORCE_TEST_PYTHON} -m unittest discover
                -s ${CMAKE_CURRENT_SOURCE_DIR}/validation_tools
                -p test_lr_force_check.py
    )
    set_tests_properties(lr_force_collector PROPERTIES TIMEOUT 30)
    add_test(
        NAME lr_energy_collector
        COMMAND ${LR_FORCE_TEST_PYTHON} -m unittest discover
                -s ${CMAKE_CURRENT_SOURCE_DIR}/validation_tools
                -p test_props_ml.py
    )
    set_tests_properties(lr_energy_collector PROPERTIES TIMEOUT 30)
endif()
