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

endif()
