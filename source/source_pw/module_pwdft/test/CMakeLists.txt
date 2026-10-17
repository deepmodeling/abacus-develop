abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_PW_td_observables
  LIBS parameter base device container psi planewave nao_smoothgen symmetry
  SOURCES
    test_td_obs_pw.cpp test_td_spin_pw.cpp test_td_vel_pw.cpp
    ../td_spin_pw.cpp ../td_vel_pw.cpp ../td_pw.cpp
    ../op_pw_vel.cpp ../projector_gradient.cpp ../velocity_workspace.cpp ../nonlocal_workspace.cpp
    ../vnl_pw.cpp ../vnl_pw_init_vnl.cpp ../vnl_pw_getvnl.cpp ../vnl_pw_deeq.cpp
    ../vnl_pw_grad.cpp ../vnl_pw_table.cpp ../vnl_pw_qrad.cpp
    ../stru_fac.cpp ../stru_fac_k.cpp ../soc.cpp ../nonlocal_maths.cpp
    ../../module_proj/onsite_proj.cpp ../../module_proj/onsite_proj_init.cpp
    ../../module_proj/onsite_proj_overlap.cpp ../../module_proj/onsite_proj_tools.cpp
    ../../module_proj/onsite_proj_tools_becp.cpp ../../module_proj/radial_proj.cpp
    ../../module_proj/onsite_proj_print.cpp
    ../../../source_basis/module_ao/orb_gaunt_table.cpp
    ../../../source_cell/klist.cpp ../../../source_cell/reciprocal_grid.cpp
    ../../../source_cell/klist_io.cpp
    ../../../source_cell/unitcell.cpp ../../../source_cell/unitcell_index.cpp
    ../../../source_cell/atom_spec.cpp ../../../source_cell/atom_pseudo.cpp ../../../source_cell/pseudo.cpp
    ../../../source_cell/magnetism.cpp ../../../source_cell/sep.cpp ../../../source_cell/sep_cell.cpp
    ../../../source_cell/read_pp.cpp ../../../source_cell/read_pp_complete.cpp
    ../../../source_cell/read_pp_upf100.cpp ../../../source_cell/read_pp_upf201.cpp
    ../../../source_cell/read_pp_vwr.cpp ../../../source_cell/read_pp_blps.cpp
    ../../../source_estate/elecstate.cpp ../../../source_estate/module_charge/charge.cpp
    ../../../source_estate/module_pot/potential_new.cpp
    ../../../source_hamilt/module_xc/xc_functional.cpp
)
target_compile_definitions(MODULE_PW_td_observables PRIVATE TD_TEST_PP_DIR="${PROJECT_SOURCE_DIR}/tests/PP_ORB")
# Keep real PW initialization and observable implementations without linking the SCF driver.
if(MSVC)
  target_compile_options(MODULE_PW_td_observables PRIVATE /Gy /Gw)
  target_link_options(MODULE_PW_td_observables PRIVATE /OPT:REF)
else()
  target_compile_options(MODULE_PW_td_observables PRIVATE -ffunction-sections -fdata-sections)
  if(APPLE)
    target_link_options(MODULE_PW_td_observables PRIVATE -Wl,-dead_strip)
  else()
    target_link_options(MODULE_PW_td_observables PRIVATE -Wl,--gc-sections)
  endif()
endif()

addtest(
  TARGET
  MODULE_PW_pwdft_soc
  LIBS
  parameter
  SOURCES
  soc_test.cpp
  ../soc.cpp
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
  TARGET
  MODULE_PW_dftu_base_test
  LIBS
  parameter
  base
  device
  SOURCES
  dftu_base_test.cpp
  ../dftu_base.cpp
  ../dftu_base_io.cpp
  ../dftu_base_io_read.cpp
  ../dftu_base_io_write.cpp
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
  TARGET
  MODULE_PW_test_dftu_pw_tools
  LIBS
  parameter
  base
  device
  SOURCES
  test_dftu_pw_tools.cpp
  ../dftu_pw_tools.cpp)

addtest(
  TARGET
  MODULE_PW_structure_factor_test
  LIBS
  parameter
  base
  device
  planewave
  SOURCES
  structure_factor_test.cpp
  ../stru_fac.cpp
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
  TARGET
  MODULE_PW_pw_hybrid_source
  LIBS
  base
  device
  SOURCES
  test_pw_hybrid_source.cpp
  ../pw_hybrid_source.cpp)
