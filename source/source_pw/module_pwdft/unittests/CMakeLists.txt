abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)


addtest(
  TARGET MODULE_PW_stress_ewa
  LIBS parameter base device planewave cell_info
  SOURCES test_stress_ewa.cpp
          ../stress_ewa.cpp
          ../../../source_cell/stru_fac.cpp
          ../../../source_cell/magnetism.cpp
          ../../../source_hamilt/module_ewald/h_ewald_pw.cpp
          ../../../source_hamilt/module_ewald/dnrm2.cpp)
# Exercise the SF shape contract in the test and its copy of stress_ewa.cpp,
# including when the surrounding ABACUS build defines NDEBUG.
target_compile_options(MODULE_PW_stress_ewa PRIVATE -UNDEBUG)
set_tests_properties(MODULE_PW_stress_ewa PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
add_test(
  NAME MODULE_PW_stress_ewa_mpi
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:MODULE_PW_stress_ewa> ${MPIEXEC_POSTFLAGS})
set_tests_properties(MODULE_PW_stress_ewa_mpi PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=1" PROCESSORS 2)

addtest(
  TARGET MODULE_PW_pwdft_soc
  LIBS parameter
  SOURCES test_soc.cpp
          ../../../source_cell/soc.cpp
          ../../../source_base/global_variable.cpp
          ../../../source_base/global_function.cpp
          ../../../source_base/global_file.cpp
          ../../../source_base/tool_quit.cpp
          ../../../source_base/global_variable.cpp
          ../../../source_base/global_file.cpp
          ../../../source_base/memory_recorder.cpp
          ../../../source_base/timer.cpp
          ../../../source_base/module_external/blas_connector_base.cpp
          ../../../source_base/module_external/blas_connector_vector.cpp
          ../../../source_base/module_external/blas_connector_matrix.cpp
          ../../../source_base/parallel_global.cpp
          ../../../source_base/parallel_comm.cpp
          ../../../source_base/parallel_common.cpp
          ../../../source_base/parallel_reduce.cpp)

addtest(
  TARGET MODULE_PW_dftu_base_test
  LIBS parameter base device
  SOURCES test_dftu_base.cpp
          ../dftu_base.cpp
          ../dftu_base_io.cpp
          ../dftu_base_read.cpp
          ../dftu_base_write.cpp
          ../dftu_pw_tools.cpp
          ../yukawa_screening.cpp
          ../../../source_estate/occ_matrix.cpp
          ../../../source_estate/occ_mixer.cpp
          ../../../source_cell/unitcell.cpp
          ../../../source_cell/unitcell_index.cpp
          ../../../source_cell/unitcell_stats.cpp
          ../../../source_cell/unitcell_setup.cpp
          ../../../source_cell/atom_spec.cpp
          ../../../source_cell/atom_pseudo.cpp
          ../../../source_cell/pseudo.cpp
          ../../../source_cell/magnetism.cpp
          ../../../source_cell/sep.cpp
          ../../../source_cell/sep_cell.cpp
          ../../../source_cell/read_atom_species.cpp
          ../../../source_cell/read_atoms.cpp
          ../../../source_cell/read_atoms_helper.cpp
          ../../../source_cell/cell_tools.cpp
          ../../../source_cell/print_cell.cpp
          ../../../source_cell/read_orb.cpp
          ../../../source_cell/read_stru.cpp
          ../../../source_cell/bcast_cell.cpp)

addtest(
  TARGET MODULE_PW_test_dftu_pw_tools
  LIBS parameter base device
  SOURCES test_dftu_pw_tools.cpp ../dftu_pw_tools.cpp)

addtest(
  TARGET MODULE_PW_structure_factor_test
  LIBS parameter base device planewave
  SOURCES test_stru_fac.cpp
          ../../../source_cell/stru_fac.cpp
          ../../../source_cell/unitcell.cpp
          ../../../source_cell/unitcell_index.cpp
          ../../../source_cell/unitcell_stats.cpp
          ../../../source_cell/unitcell_setup.cpp
          ../../../source_cell/update_cell.cpp
          ../../../source_cell/bcast_cell.cpp
          ../../../source_cell/print_cell.cpp
          ../../../source_cell/atom_spec.cpp
          ../../../source_cell/atom_pseudo.cpp
          ../../../source_cell/pseudo.cpp
          ../../../source_cell/read_stru.cpp
          ../../../source_cell/read_atom_species.cpp
          ../../../source_cell/read_atoms.cpp
          ../../../source_cell/read_atoms_helper.cpp
          ../../../source_cell/cell_tools.cpp
          ../../../source_cell/read_pp.cpp
          ../../../source_cell/read_pp_complete.cpp
          ../../../source_cell/read_pp_upf100.cpp
          ../../../source_cell/read_pp_upf201.cpp
          ../../../source_cell/read_pp_vwr.cpp
          ../../../source_cell/read_pp_blps.cpp
          ../../../source_cell/sep.cpp
          ../../../source_cell/sep_cell.cpp
          ../../../source_cell/read_pp_ucell.cpp
          ../../../source_cell/cal_wfc.cpp
          ../../../source_cell/cal_nelec_nband.cpp
          ../../../source_cell/read_orb.cpp)

addtest(
  TARGET MODULE_PW_pw_hybrid_source
  LIBS base device
  SOURCES test_pw_hybrid_source.cpp ../pw_hybrid_source.cpp)
