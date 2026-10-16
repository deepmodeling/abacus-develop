abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

if (ENABLE_MPI)

AddTest(
  TARGET MODULE_CHARGE_rhog_io
  LIBS parameter base device planewave
  SOURCES test_chg_rhog_io.cpp ../chg_rhog_io.cpp ../../../source_basis/module_pw/test/test_tool.cpp
  # This test drives PW_Basis::initmpi and read/write_rhog's MPI collectives,
  # so it must keep __MPI even though this directory disables it. Its main()
  # calls MPI_Init via test_tool.cpp's setupmpi().
  KEEP_FEATURE_DEFINITIONS __MPI
)

AddTest(
  TARGET MODULE_CHARGE_charge
  LIBS parameter planewave_serial base device cell_info
  SOURCES test_charge.cpp ../charge.cpp ../chg_atomic.cpp ../chg_atomic_inner.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_extra
  LIBS parameter base device cell_info
  SOURCES test_chg_extra.cpp ../chg_extra.cpp ../../../source_io/module_output/read_cube.cpp ../../../source_io/module_output/write_cube.cpp
 ../../../source_base/module_fft/fft_bundle.cpp ../../../source_base/module_fft/fft_cpu.cpp
)

AddTest(
  TARGET MODULE_CHARGE_tools
  LIBS parameter planewave_serial base device cell_info
  SOURCES test_chg_tools.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_uspp
  LIBS base device
  SOURCES test_chg_uspp.cpp ../chg_uspp.cpp
)

AddTest(
  TARGET MODULE_CHARGE_dmr
  LIBS base device
  SOURCES test_chg_dmr.cpp ../chg_dmr.cpp
)

AddTest(
  TARGET MODULE_CHARGE_precond
  LIBS base device planewave_serial
  SOURCES test_chg_precond.cpp ../chg_precond.cpp
)

AddTest(
  TARGET MODULE_CHARGE_drho
  LIBS parameter planewave_serial base device cell_info
  SOURCES test_chg_drho.cpp ../chg_drho.cpp ../chg_drho_inner.cpp ../charge.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_drho_inner
  LIBS base device planewave_serial
  SOURCES test_chg_drho_inner.cpp ../chg_drho_inner.cpp
)

AddTest(
  TARGET MODULE_CHARGE_mix_rho
  LIBS parameter base psi device planewave_serial cell_info
  SOURCES test_chg_mix_rho.cpp
  ../chg_mix.cpp ../chg_mix_rho.cpp
  ../chg_dmr.cpp
  ../chg_drho.cpp
  ../chg_drho_inner.cpp
  ../chg_precond.cpp ../chg_tau.cpp
  ../chg_uspp.cpp
  ../charge.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_mix
  LIBS parameter base  psi  device planewave_serial cell_info
  SOURCES test_chg_mix.cpp
  ../chg_mix.cpp ../chg_mix_rho.cpp
  ../chg_dmr.cpp
  ../chg_drho.cpp
  ../chg_drho_inner.cpp
  ../chg_precond.cpp ../chg_tau.cpp
  ../chg_uspp.cpp
)

AddTest(
  TARGET MODULE_CHARGE_symm
  LIBS parameter base device cell_info planewave_serial symmetry
  SOURCES test_chg_symm.cpp ../chg_symm.cpp ../chg_symm_detail.cpp ../charge.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_symm_detail
  LIBS base device cell_info planewave_serial symmetry
  SOURCES test_chg_symm_detail.cpp ../chg_symm_detail.cpp
)

AddTest(
  TARGET MODULE_CHARGE_atomic_inner
  LIBS parameter base device cell_info planewave_serial
  SOURCES test_chg_atomic_inner.cpp ../chg_atomic_inner.cpp ../charge.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_atomic
  LIBS parameter base device cell_info planewave_serial
  SOURCES test_chg_atomic.cpp ../chg_atomic.cpp ../chg_atomic_inner.cpp ../charge.cpp ../chg_tools.cpp
)

AddTest(
  TARGET MODULE_CHARGE_tau
  LIBS parameter base psi device planewave_serial cell_info
  SOURCES test_chg_tau.cpp ../chg_tau.cpp ../chg_uspp.cpp ../charge.cpp ../chg_tools.cpp
)

# MODULE_CHARGE_routine is temporarily disabled: chg_routine.cpp couples to
# heavy external symbols (Plus_U_Base from source_pw/module_pwdft, elecstate
# ::update_pot / ElecState::cal_energies) whose transitive dependencies are
# hard to satisfy in this test; chg_mix_rho.cpp also calls
# module_charge::detail::mix_tau_recip so chg_tau.cpp must be linked too.
# TODO: untangle this deep coupling later before re-enabling.
# AddTest(
#   TARGET MODULE_CHARGE_routine
#   LIBS parameter base psi device planewave_serial cell_info
#   SOURCES test_chg_routine.cpp ../chg_routine.cpp ../chg_dmr.cpp ../charge.cpp ../chg_tools.cpp
#   ../chg_mix.cpp ../chg_mix_rho.cpp ../chg_drho.cpp ../chg_drho_inner.cpp ../chg_precond.cpp ../chg_uspp.cpp
# )

# MODULE_CHARGE_init is temporarily disabled: chg_init.cpp pulls in symbols
# from source_io (ModuleIO::read_vdata_palgrid, ModuleIO::read_wf2rho_pw,
# Restart::read_file2) whose dependency chains are complex (read_wf2rho_pw
# alone drags in chg_parallel/chg_symm/elecstate_op), and the test file also
# stubs GlobalC::restart which conflicts with linking the real restart.cpp.
# TODO: resolve this dependency issue later before re-enabling.
# AddTest(
#   TARGET MODULE_CHARGE_init
#   LIBS parameter base device cell_info planewave_serial
#   SOURCES test_chg_init.cpp ../chg_init.cpp ../chg_atomic.cpp ../chg_atomic_inner.cpp ../chg_rhog_io.cpp ../charge.cpp ../chg_tools.cpp
# )

AddTest(
  TARGET MODULE_CHARGE_parallel
  LIBS parameter psi base device planewave
  SOURCES test_chg_parallel.cpp ../chg_parallel.cpp
  # Real MPI test: its main() calls MPI_Init unconditionally, so it must keep
  # __MPI despite this directory disabling it.
  KEEP_FEATURE_DEFINITIONS __MPI
)

add_test(NAME MODULE_CHARGE_parallel_4np
      COMMAND mpirun -np 4 ./MODULE_CHARGE_parallel;
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

endif()
