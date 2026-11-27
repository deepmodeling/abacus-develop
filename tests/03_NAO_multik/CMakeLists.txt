enable_testing()

find_program(BASH bash)
if(ENABLE_ASAN)
  add_test(
    NAME 03_NAO_multik_test_with_asan
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 2 -s true
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
else()
  add_test(
    NAME 03_NAO_multik
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
endif()

# The consumer check requires the external Wannier90 and postw90 executables.
find_program(WANNIER90_EXECUTABLE wannier90.x)
find_program(POSTW90_EXECUTABLE postw90.x)
if(ENABLE_LCAO
   AND ENABLE_MPI
   AND WANNIER90_EXECUTABLE
   AND POSTW90_EXECUTABLE)
  find_package(Python3 QUIET COMPONENTS Interpreter NumPy)
  if(Python3_NumPy_FOUND AND MPIEXEC_EXECUTABLE)
    add_test(
      NAME 03_NAO_wannier_spn
      COMMAND
        ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/interfaces/Wannier90_interface/tests/check_wannier_spn.py
        ${WANNIER90_EXECUTABLE} ${POSTW90_EXECUTABLE} ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${ABACUS_BIN_PATH}
        ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      03_NAO_wannier_spn PROPERTIES PROCESSORS 4 TIMEOUT 300 ENVIRONMENT
                                    "OMP_NUM_THREADS=1")
  endif()
endif()

if(ENABLE_CNPY)
  add_test(
    NAME 03_NAO_multik_npz
    COMMAND ${BASH} ../integrate/Autotest.sh -a ${ABACUS_BIN_PATH} -n 4 -j auto
            -f CASES_CNPY.txt
    WORKING_DIRECTORY ${ABACUS_TEST_DIR}/03_NAO_multik)
endif()
