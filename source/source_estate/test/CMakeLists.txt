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
  TARGET MODULE_ESTATE_Elecstate_Op_UTs
  LIBS parameter  psi base device
  SOURCES ../kernels/test/elecstate_op_test.cpp
)

if(USE_CUDA)
  target_compile_definitions(MODULE_ESTATE_Elecstate_Op_UTs PRIVATE __UT_USE_CUDA)
elseif(USE_ROCM)
  target_compile_definitions(MODULE_ESTATE_Elecstate_Op_UTs PRIVATE __UT_USE_ROCM)
endif()

AddTest(
  TARGET MODULE_ESTATE_elecstate_occupy
  LIBS parameter  base device
  SOURCES elecstate_occupy_test.cpp ../occupy.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_fp_energy
  LIBS parameter  base device
  SOURCES elecstate_fp_energy_test.cpp ../fp_energy.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_print
  LIBS parameter  base device symmetry
  SOURCES elecstate_print_test.cpp ../elecstate_print.cpp ../occupy.cpp
  ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_base
  LIBS parameter  base device symmetry
  SOURCES elecstate_base_test.cpp ../elecstate.cpp  ../elecstate_tools.cpp ../occupy.cpp ../../source_psi/psi.cpp
  ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
  ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_pw
  LIBS parameter  planewave_serial base device symmetry
  SOURCES elecstate_pw_test.cpp
    ../elecstate_pw.cpp
    ../elecstate_pw_cal_tau.cpp
    ../elecstate.cpp
    ../occupy.cpp
    ../module_charge/charge_mpi.cpp
    ../../source_lcao/module_deltaspin/spin_constrain.cpp
    ../../source_lcao/module_deltaspin/deltaspin_state.cpp
    ../../source_psi/psi.cpp
    ../../source_base/module_device/memory_op.cpp
    ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_elecstate_energy
  LIBS parameter  base device planewave_serial symmetry
  SOURCES elecstate_energy_test.cpp
    ../elecstate_energy.cpp
    ../fp_energy.cpp
    ../makov_payne.cpp
    ../module_pot/h_hartree_pw.cpp
    ../../source_cell/klist.cpp ../../source_cell/klist_io.cpp ../../source_cell/parallel_kpoints.cpp ../../source_cell/reciprocal_grid.cpp
)

AddTest(
  TARGET MODULE_ESTATE_potentials_base
  SOURCES potentials_base_test.cpp
)

AddTest(
  TARGET MODULE_ESTATE_occ_mixer
  LIBS parameter base device cell_info
  SOURCES test_occ_mixer.cpp ../occ_mixer.cpp ../occ_matrix.cpp
)

AddTest(
  TARGET MODULE_ESTATE_potentials_new
  LIBS parameter  base device planewave_serial
  SOURCES potential_new_test.cpp td_field_test.cpp ../module_pot/potential_new.cpp
    ../module_pot/td_field.cpp ../module_pot/td_field_profiles.cpp ../module_pot/td_field_manager.cpp
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
  ../module_charge/charge_mixing.cpp ../module_charge/charge_mixing_dmr.cpp ../module_charge/mix_resid.cpp
  ../module_charge/mix_precond.cpp ../module_charge/charge_mixing_rho.cpp
  ../module_charge/charge_mixing_uspp.cpp
)

AddTest(
  TARGET MODULE_ESTATE_charge_extra
  LIBS parameter base device cell_info
  SOURCES charge_extra_test.cpp ../module_charge/charge_extra.cpp ../../source_io/module_output/read_cube.cpp ../../source_io/module_output/write_cube.cpp
 ../../source_base/module_fft/fft_bundle.cpp ../../source_base/module_fft/fft_cpu.cpp
)

AddTest(
  TARGET MODULE_ESTATE_gint_prec_ctrl
  LIBS parameter base device
  SOURCES gint_prec_ctrl_test.cpp
    ../module_charge/gint_prec_ctrl.cpp
)

AddTest(
  TARGET MODULE_ESTATE_test_rhog_io
  LIBS parameter base device planewave
  SOURCES test_rhog_io.cpp ../rhog_io.cpp ../../source_basis/module_pw/test/test_tool.cpp
  # This test drives PW_Basis::initmpi and read/write_rhog's MPI collectives,
  # so it must keep __MPI even though this directory disables it. Its main()
  # calls MPI_Init via test_tool.cpp's setupmpi().
  KEEP_FEATURE_DEFINITIONS __MPI
)

AddTest(
  TARGET MODULE_ESTATE_charge_mpi_test
  LIBS parameter psi base device planewave
  SOURCES charge_mpi_test.cpp ../module_charge/charge_mpi.cpp
  # Real MPI test: its main() calls MPI_Init unconditionally, so it must keep
  # __MPI despite this directory disabling it.
  KEEP_FEATURE_DEFINITIONS __MPI
)

add_test(NAME MODULE_ESTATE_charge_mpi_test_4np
      COMMAND mpirun -np 4 ./MODULE_ESTATE_charge_mpi_test;
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

endif()
