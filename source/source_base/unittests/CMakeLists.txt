abacus_disable_feature_definitions(__MPI)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_BASE_csr_reader
  LIBS parameter base device
  SOURCES test_csr_file_reader.cpp test_csr_reader.cpp)

AddTest(
  TARGET MODULE_BASE_filename
  LIBS parameter base device
  SOURCES test_filename.cpp
)

AddTest(
  TARGET MODULE_BASE_matrix2d
  LIBS parameter base device container
  SOURCES test_matrix2d.cpp
  KEEP_FEATURE_DEFINITIONS __MPI
)

AddTest(
  TARGET MODULE_BASE_matrix2d_linalg
  LIBS parameter base device container
  SOURCES test_matrix2d_linalg.cpp
  KEEP_FEATURE_DEFINITIONS __MPI
)

set_tests_properties(
  MODULE_BASE_matrix2d
  MODULE_BASE_matrix2d_linalg
  PROPERTIES ENVIRONMENT
             "OMP_NUM_THREADS=1;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1;HWLOC_COMPONENTS=-gl"
)

if(ENABLE_MPI)
  foreach(matrix2d_target IN ITEMS
          MODULE_BASE_matrix2d
          MODULE_BASE_matrix2d_linalg)
    get_target_property(
      matrix2d_kept_feature_definitions
      ${matrix2d_target}
      ABACUS_KEPT_FEATURE_DEFINITIONS)
    list(FIND matrix2d_kept_feature_definitions "__MPI" mpi_feature_index)
    if(mpi_feature_index EQUAL -1)
      message(
        FATAL_ERROR
        "${matrix2d_target} must keep __MPI to match the linked base library")
    endif()
  endforeach()

  add_test(
    NAME MODULE_BASE_matrix2d_2ranks
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_BASE_matrix2d>
            ${MPIEXEC_POSTFLAGS}
  )
  add_test(
    NAME MODULE_BASE_matrix2d_4ranks
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_BASE_matrix2d>
            ${MPIEXEC_POSTFLAGS}
  )
  add_test(
    NAME MODULE_BASE_matrix2d_linalg_2ranks
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_BASE_matrix2d_linalg>
            ${MPIEXEC_POSTFLAGS}
  )
  add_test(
    NAME MODULE_BASE_matrix2d_linalg_4ranks
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_BASE_matrix2d_linalg>
            ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    MODULE_BASE_matrix2d_2ranks
    MODULE_BASE_matrix2d_4ranks
    MODULE_BASE_matrix2d_linalg_2ranks
    MODULE_BASE_matrix2d_linalg_4ranks
    PROPERTIES ENVIRONMENT
               "OMP_NUM_THREADS=1;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1;HWLOC_COMPONENTS=-gl"
  )
endif()
