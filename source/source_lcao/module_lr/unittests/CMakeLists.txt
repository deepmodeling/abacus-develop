if(ENABLE_LIBXC)
  addtest(
    TARGET MODULE_LR_pot_hxc_lrtd
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_pot_hxc_lrtd.cpp
            ../potentials/pot_hxc_lrtd.cpp
            ../kernel_pots.cpp
            ../potentials/pot_grad_xc.cpp
            ../potentials/xc_kernel.cpp
            ../potentials/xc_kxc.cpp
            ../utils/lr_util.cpp
            ../../../source_hamilt/module_xc/xc_functional.cpp
            ../../../source_hamilt/module_xc/libxc_setup.cpp)

  addtest(
    TARGET MODULE_LR_xc_kernel
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_xc_kernel.cpp
            ../potentials/xc_kxc.cpp
            ../potentials/xc_kernel.cpp
            ../utils/lr_util.cpp
            ../../../source_hamilt/module_xc/xc_functional.cpp
            ../../../source_hamilt/module_xc/libxc_setup.cpp)

  addtest(
    TARGET MODULE_LR_pot_grad_xc
    LIBS base parameter device container psi Libxc::xc
    SOURCES test_pot_grad_xc.cpp
            ../potentials/pot_grad_xc.cpp
            ../potentials/xc_kernel.cpp
            ../potentials/xc_kxc.cpp
            ../utils/lr_util.cpp
            ../../../source_hamilt/module_xc/xc_functional.cpp
            ../../../source_hamilt/module_xc/libxc_setup.cpp)

  if(ENABLE_MPI)
    foreach(target MODULE_LR_pot_hxc_lrtd MODULE_LR_xc_kernel MODULE_LR_pot_grad_xc)
      target_sources(${target} PRIVATE ../../../source_base/test/mpi_test_main.cpp)
    endforeach()
    add_test(NAME MODULE_LR_xc_kernel_parallel
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:MODULE_LR_xc_kernel> ${MPIEXEC_POSTFLAGS}
        --gtest_filter=KernelXCTest.HseGroundResponseFallback*:KernelXCTest.HseGroundResponseRetry*)
    set_tests_properties(MODULE_LR_xc_kernel_parallel PROPERTIES
      PROCESSORS 4 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
  endif()
endif()
