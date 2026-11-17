addtest(
  TARGET MODULE_LR_nacv_result
  SOURCES test_nacv_result.cpp ../nacv_result.cpp ../../../source_base/matrix.cpp)
# This isolated result test does not exercise matrix BLAS operations.
target_compile_definitions(MODULE_LR_nacv_result PRIVATE __NORMAL)

addtest(
  TARGET MODULE_LR_nacv_io
  LIBS base parameter device container psi
  SOURCES test_nacv_io.cpp ../nacv_io.cpp)
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
    TARGET MODULE_LR_nacv_cg
    LIBS base parameter container device psi MPI::MPI_CXX
    SOURCES test_nacv_cg.cpp ../nacv_cg.cpp ../../../source_base/test/mpi_test_main.cpp)
  add_test(
    NAME MODULE_LR_nacv_cg_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_LR_nacv_cg> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_LR_nacv_cg_mpi4 PROPERTIES PROCESSORS 4 TIMEOUT 120)
  set_tests_properties(MODULE_LR_nacv_cg MODULE_LR_nacv_cg_mpi4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

addtest(
  TARGET MODULE_LR_nacv_xc_names
  SOURCES test_nacv_xc_names.cpp)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_LR_nacv_factorization
    LIBS base parameter ${math_libs} container device psi MPI::MPI_CXX
    SOURCES test_nacv_factorization.cpp ../nacv_factorization.cpp ../utils/lr_util.cpp
            ../../../source_base/test/mpi_test_main.cpp)
  add_test(
    NAME MODULE_LR_nacv_factorization_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_LR_nacv_factorization> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_LR_nacv_factorization_mpi4 PROPERTIES PROCESSORS 4 TIMEOUT 120)
  set_tests_properties(MODULE_LR_nacv_factorization MODULE_LR_nacv_factorization_mpi4
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
