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
if(ENABLE_LIBXC AND ENABLE_MPI)
    add_test(
        NAME PW_HYBRID_NSCF
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_hybrid_nscf.py
                ${ABACUS_BIN_PATH} --mpi-ranks 2
    )
    set_tests_properties(PW_HYBRID_NSCF PROPERTIES
        LABELS "PW;EXX;NSCF" PROCESSORS 2 TIMEOUT 300
    )
endif()

if(ENABLE_LIBXC AND ENABLE_MPI AND USE_CUDA)
    add_test(
        NAME PW_HYBRID_NSCF_GPU
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_hybrid_nscf.py
                ${ABACUS_BIN_PATH} --gpu
    )
    set_tests_properties(PW_HYBRID_NSCF_GPU PROPERTIES
        LABELS "PW;EXX;NSCF;GPU" PROCESSORS 1 TIMEOUT 300 RUN_SERIAL TRUE
    )
endif()
