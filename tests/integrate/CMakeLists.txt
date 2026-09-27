find_program(BASH bash)
find_package(Python3 REQUIRED)

set(_dftb_native_mpi_args)
if(ENABLE_MPI)
    list(APPEND _dftb_native_mpi_args
        --mpi-exec=${MPIEXEC_EXECUTABLE}
        --mpi-num-procs-flag=${MPIEXEC_NUMPROC_FLAG}
    )
    foreach(_mpi_flag IN LISTS MPIEXEC_PREFLAGS)
        list(APPEND _dftb_native_mpi_args --mpi-preflag=${_mpi_flag})
    endforeach()
    foreach(_mpi_flag IN LISTS MPIEXEC_POSTFLAGS)
        list(APPEND _dftb_native_mpi_args --mpi-postflag=${_mpi_flag})
    endforeach()
endif()

add_test(
    NAME dftb_native_c2n_reference
    COMMAND ${Python3_EXECUTABLE}
        ${ABACUS_TEST_DIR}/03_DFTB_native_C2N/run_regression.py
        --abacus ${ABACUS_BIN_PATH}
        --case-dir ${ABACUS_TEST_DIR}/03_DFTB_native_C2N
        ${_dftb_native_mpi_args}
)
set_tests_properties(dftb_native_c2n_reference PROPERTIES
    TIMEOUT 900
    LABELS "integration;dftb"
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
