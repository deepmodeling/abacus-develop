if(ENABLE_MPI)
  set(test_cvcx_main ../../../../source_base/test/mpi_test_main.cpp)
endif()
abacus_disable_feature_definitions(__LIBXC)
addtest(
  TARGET MODULE_LR_ao_to_mo_test
  LIBS parameter base container device psi
  SOURCES ao_to_mo_test.cpp ../../utils/lr_util.cpp ../ao_to_mo_parallel.cpp
          ../ao_to_mo_serial.cpp)

addtest(
  TARGET MODULE_LR_CVCX_test
  LIBS base parameter ${math_libs} container device psi
  SOURCES test_cvcx.cpp ${test_cvcx_main} ../../utils/lr_util.cpp
          ../cvcx_par.cpp ../cvcx_serial.cpp)
if(ENABLE_MPI)
  addtest(
    TARGET MODULE_LR_empty_blocks
    LIBS parameter base container device psi
    SOURCES test_ao_to_mo_parallel.cpp
            ${test_cvcx_main}
            ../../utils/lr_util.cpp
            ../ao_to_mo_parallel.cpp
            ../cvcx_par.cpp
            ../cvcx_serial.cpp
            ../../dm_trans/dm_trans_parallel.cpp
            ../../dm_trans/dm_trans_serial.cpp
            ../../dm_trans/dm_diff_par.cpp
            ../../dm_trans/dm_diff_ser.cpp)
  add_test(NAME MODULE_LR_empty_blocks_mpi
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:MODULE_LR_empty_blocks>)
  set_tests_properties(
    MODULE_LR_empty_blocks_mpi
    PROPERTIES PROCESSORS 4 TIMEOUT 60 ENVIRONMENT
               "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
