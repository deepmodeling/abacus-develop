abacus_disable_feature_definitions(__MPI __LCAO)
abacus_add_local_feature_definitions(__NORMAL)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_MD_dump_mpi
    KEEP_FEATURE_DEFINITIONS __MPI
    LIBS parameter psi device
    SOURCES test_md_func.cpp ../../source_base/parallel_common.cpp
            ../../source_cell/magnetism.cpp ${md_test_depend_files})
  add_test(
    NAME MODULE_MD_dump_mpi_np4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_MD_dump_mpi> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(
    MODULE_MD_dump_mpi MODULE_MD_dump_mpi_np4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;I_MPI_PIN=0" TIMEOUT 60)
endif()
