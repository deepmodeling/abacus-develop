addtest(
  TARGET MODULE_LR_nacv_result
  SOURCES test_nacv_result.cpp ../nacv_result.cpp ../../../source_base/matrix.cpp)
# This isolated result test does not exercise matrix BLAS operations.
target_compile_definitions(MODULE_LR_nacv_result PRIVATE __NORMAL)

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
    TARGET MODULE_LR_zeq_cg_precision
    LIBS base parameter container device psi MPI::MPI_CXX
    SOURCES test_zeq_cg.cpp ../zeq_cg.cpp ../../../source_base/test/mpi_test_main.cpp)
  add_test(
    NAME MODULE_LR_zeq_cg_precision_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_LR_zeq_cg_precision> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_LR_zeq_cg_precision_mpi4 PROPERTIES PROCESSORS 4 TIMEOUT 120)
  set_tests_properties(MODULE_LR_zeq_cg_precision MODULE_LR_zeq_cg_precision_mpi4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
