abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

if (ENABLE_MPI)

install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

AddTest(
  TARGET MODULE_ESTATE_Elecstate_Op_UTs
  LIBS parameter  psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_occupy
  LIBS parameter  base device
  SOURCES elecstate_occupy_test.cpp ../occupy.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_magnetism
  LIBS parameter  base device
  SOURCES elecstate_magnetism_test.cpp ../../source_cell/magnetism.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_fp_energy
  LIBS parameter  base device
  SOURCES elecstate_fp_energy_test.cpp ../fp_energy.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_print
  LIBS parameter  base device
  SOURCES elecstate_print_test.cpp ../elecstate_print.cpp ../occupy.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_base
  LIBS parameter  base device
  SOURCES elecstate_base_test.cpp ../elecstate.cpp  ../elecstate_tools.cpp ../occupy.cpp ../../source_psi/psi.cpp
  ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_pw
  LIBS parameter  planewave_serial base device
  SOURCES elecstate_pw_test.cpp
    ../elecstate_pw.cpp
    ../elecstate_pw_cal_tau.cpp
    ../elecstate.cpp
    ../occupy.cpp
    ../module_charge/charge_mpi.cpp
    ../../source_psi/psi.cpp
    ../../source_base/module_device/memory_op.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_energy
  LIBS parameter  base device planewave_serial
  SOURCES elecstate_energy_test.cpp
    ../elecstate_energy.cpp
    ../fp_energy.cpp
    ../makov_payne.cpp
    ../module_pot/H_Hartree_pw.cpp
)

AddTest(
  TARGET MODULE_ESTATE_potentials_base
  SOURCES potentials_base_test.cpp
)

AddTest(
  TARGET MODULE_ESTATE_potentials_new
  LIBS parameter  base device planewave_serial
  SOURCES potential_new_test.cpp ../module_pot/potential_new.cpp
)

AddTest(
  TARGET MODULE_ESTATE_charge_test
  LIBS parameter  planewave_serial base device cell_info
  SOURCES charge_test.cpp ../module_charge/charge.cpp

)

AddTest(
  TARGET MODULE_ESTATE_charge_mixing
  LIBS parameter base  psi  device planewave_serial cell_info
  SOURCES charge_mixing_test.cpp  
  ../module_charge/charge_mixing.cpp ../module_charge/charge_mixing_dmr.cpp ../module_charge/charge_mixing_residual.cpp
  ../module_charge/charge_mixing_preconditioner.cpp ../module_charge/charge_mixing_rho.cpp
  ../module_charge/charge_mixing_uspp.cpp
)

AddTest(
  TARGET MODULE_ESTATE_charge_extra
  LIBS parameter base device cell_info
  SOURCES charge_extra_test.cpp ../module_charge/charge_extra.cpp ../../source_io/module_output/read_cube.cpp ../../source_io/module_output/write_cube.cpp
 ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
)

AddTest(
  TARGET MODULE_ESTATE_gint_precision_controller
  LIBS parameter base device
  SOURCES gint_precision_controller_test.cpp
    ../module_charge/gint_precision_controller.cpp
)

endif()
