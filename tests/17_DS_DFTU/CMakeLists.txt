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
        COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
        WORKING_DIRECTORY ${ABACUS_TEST_DIR}/17_DS_DFTU
    )
endif()

if(ENABLE_LCAO AND ENABLE_LIBXC AND ENABLE_MPI)
    find_package(Python3 REQUIRED)
    add_test(
        NAME lcao_deltaspin_force
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../integrate/test_deltaspin_force.py
            --abacus ${ABACUS_BIN_PATH}
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/65_LCAO_DS_S4_SO_CF
            --device cpu
            --launcher "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4"
    )
    set_tests_properties(lcao_deltaspin_force PROPERTIES
        LABELS "integration;lcao;deltaspin"
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1"
        PROCESSORS 4
        TIMEOUT 900
        RUN_SERIAL TRUE
    )
endif()
