addtest(
  TARGET MODULE_LR_nacv_tools
  LIBS base parameter ${math_libs} device container psi
  SOURCES test_nacv_tools.cpp test_lr_util_xc.cpp ../nacv_tools.cpp ../utils/lr_util.cpp)
# State serialization must also execute with assertions disabled.
target_compile_definitions(MODULE_LR_nacv_tools PRIVATE NDEBUG)

if(ENABLE_MPI)
  target_sources(MODULE_LR_nacv_tools PRIVATE ../../../source_base/test/mpi_test_main.cpp)
  target_link_libraries(MODULE_LR_nacv_tools PRIVATE MPI::MPI_CXX)
  add_test(
    NAME MODULE_LR_nacv_tools_mpi4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:MODULE_LR_nacv_tools>
            --gtest_filter=ZeqCgPrecision.*:NacvFactorization.*:NacvIoMpi.* ${MPIEXEC_POSTFLAGS})
  set_tests_properties(MODULE_LR_nacv_tools_mpi4 PROPERTIES PROCESSORS 4 TIMEOUT 120)
  # Both invocations of NacvIoMpi use the same scratch filenames.
  set_tests_properties(MODULE_LR_nacv_tools MODULE_LR_nacv_tools_mpi4
    PROPERTIES RESOURCE_LOCK nacv_tools_scratch
               ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

if(ENABLE_LIBXC)
  addtest(
    TARGET MODULE_LR_xc_kernel
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_xc_kernel.cpp ../potentials/xc_kernel.cpp ../potentials/xc_kxc.cpp
            ../utils/lr_util.cpp ../../../source_hamilt/module_xc/xc_functional.cpp
            ../../../source_hamilt/module_xc/libxc_setup.cpp)
endif()
