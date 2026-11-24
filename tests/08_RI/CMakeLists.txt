enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
  add_test(
    NAME 08_RI_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI)
else()
  add_test(
    NAME 08_RI
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI)

  # These compact cases validate only producer artifacts consumed by LibRPA.
  # They leave the established EXX and BSE regression corpus unchanged.
  add_test(
    NAME 08_RI_librpa_producer
    COMMAND ${BASH} ../integrate/validation_tools/run_librpa_producer.sh
            ${ABACUS_BIN_PATH} ${ABACUS_TEST_DIR}/08_RI
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/08_RI)
  # Keep launcher flags in one argv item because they start with a dash and
  # the Python helper's argparse would otherwise parse their values.
  set(MPI_LAUNCHER_FLAGS)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND MPI_LAUNCHER_FLAGS "--mpi-preflag=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND MPI_LAUNCHER_FLAGS "--mpi-postflag=${flag}")
  endforeach()
  add_test(
    NAME 08_RI_librpa_wavefunction_checks
    COMMAND
      ${BASH}
      ${ABACUS_TEST_DIR}/integrate/validation_tools/run_librpa_wavefunction_checks.sh
      ${ABACUS_TEST_DIR}/integrate/validation_tools)
  add_test(
    NAME 08_RI_librpa_producer_mpi
    COMMAND
      ${BASH} ${ABACUS_TEST_DIR}/integrate/validation_tools/run_librpa_mpi.sh
      --abacus ${ABACUS_BIN_PATH} --cases ${ABACUS_TEST_DIR}/08_RI --mpirun
      ${MPIEXEC_EXECUTABLE} --mpi-np-flag=${MPIEXEC_NUMPROC_FLAG}
      ${MPI_LAUNCHER_FLAGS})
  set_tests_properties(08_RI_librpa_producer_mpi PROPERTIES PROCESSORS 2
                                                            TIMEOUT 1200)
  add_test(
    NAME 08_RI_librpa_symmetry
    COMMAND
      ${BASH} ${ABACUS_TEST_DIR}/integrate/validation_tools/run_librpa_sym.sh
      --abacus ${ABACUS_BIN_PATH} --cases ${ABACUS_TEST_DIR}/08_RI --mpirun
      ${MPIEXEC_EXECUTABLE} --mpi-np-flag=${MPIEXEC_NUMPROC_FLAG}
      ${MPI_LAUNCHER_FLAGS})
  set_tests_properties(08_RI_librpa_symmetry PROPERTIES PROCESSORS 2 TIMEOUT
                                                        1800)

  # This two-stage producer/consumer check needs the same SCF artifacts in
  # separate LR runs, which the single-invocation Autotest cases cannot express.
  find_program(LR_ABSORPTION_PYTHON python3)
  if(ENABLE_MPI AND ENABLE_LCAO AND ENABLE_LIBRI AND ENABLE_LIBXC AND LR_ABSORPTION_PYTHON)
    add_test(
      NAME 08_RI_lr_absorption_workflow
      COMMAND
        ${LR_ABSORPTION_PYTHON}
        ${ABACUS_TEST_DIR}/integrate/validation_tools/check_lr_absorption_workflow.py
        --abacus ${ABACUS_BIN_PATH}
        --case ${ABACUS_TEST_DIR}/08_RI/lr_tddft_velocity_commutator
        --mpiexec ${MPIEXEC_EXECUTABLE} --mpi-np-flag=${MPIEXEC_NUMPROC_FLAG}
        ${MPI_LAUNCHER_FLAGS})
    set_tests_properties(08_RI_lr_absorption_workflow PROPERTIES PROCESSORS 4 TIMEOUT 180)
  endif()
endif()
