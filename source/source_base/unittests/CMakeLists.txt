AddTest(
  TARGET MODULE_BASE_matrix2d
  LIBS parameter base device container
  SOURCES test_matrix2d.cpp
)

AddTest(
  TARGET MODULE_BASE_matrix2d_linalg
  LIBS parameter base device container
  SOURCES test_matrix2d_linalg.cpp
)

set_tests_properties(
  MODULE_BASE_matrix2d
  MODULE_BASE_matrix2d_linalg
  PROPERTIES ENVIRONMENT
             "OMP_NUM_THREADS=1;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1;HWLOC_COMPONENTS=-gl"
)

if(ENABLE_MPI)
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
