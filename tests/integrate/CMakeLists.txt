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

if(ENABLE_LCAO AND ENABLE_LIBXC AND ENABLE_MPI)
    add_test(
        NAME lcao_deltaspin_force
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_deltaspin_force.py
            --abacus ${ABACUS_BIN_PATH}
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/deltaspin_force
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
