addtest(
  TARGET MODULE_LR_nacv_result
  LIBS base parameter ${math_libs} device container psi
  SOURCES test_nacv_tools_result.cpp ../nacv_tools.cpp ../utils/lr_util.cpp)

addtest(
  TARGET MODULE_LR_nacv_io
  LIBS base parameter ${math_libs} device container psi
  SOURCES test_nacv_tools_io.cpp ../nacv_tools.cpp ../utils/lr_util.cpp)
# State export must run even when assertions are disabled.
target_compile_definitions(MODULE_LR_nacv_io PRIVATE NDEBUG)

if(ENABLE_LIBXC)
  addtest(
    TARGET MODULE_LR_xc_kernel
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_xc_kernel.cpp ../potentials/xc_kernel.cpp ../potentials/xc_kxc.cpp
            ../utils/lr_util.cpp ../../../source_hamilt/module_xc/xc_functional.cpp
            ../../../source_hamilt/module_xc/libxc_setup.cpp)
endif()

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_LR_nacv_solver
    LIBS base parameter ${math_libs} container device psi MPI::MPI_CXX
    SOURCES test_nacv_tools_solver.cpp ../nacv_tools.cpp ../utils/lr_util.cpp
            ../../../source_base/test/mpi_test_main.cpp)
  add_test(
    NAME MODULE_LR_nacv_solver_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_LR_nacv_solver> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_LR_nacv_solver_mpi4 PROPERTIES PROCESSORS 4 TIMEOUT 120)
  set_tests_properties(MODULE_LR_nacv_solver MODULE_LR_nacv_solver_mpi4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

addtest(
  TARGET MODULE_LR_nacv_xc_names
  SOURCES test_nacv_tools_xc_names.cpp)
