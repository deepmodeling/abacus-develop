abacus_disable_feature_definitions(__MPI)
abacus_disable_feature_definitions(__EXX)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__UT_USE_CUDA)
abacus_disable_feature_definitions(__UT_USE_ROCM)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(_OPENMP)

if (ENABLE_MPI)

# Copy at configure time so a plain `make` + `ctest` run finds the data.
# install() only runs during `cmake --install`, which local test runs skip.
file(COPY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

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
