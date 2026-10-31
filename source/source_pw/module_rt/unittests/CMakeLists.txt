abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

set(LINALG_DIR ${PROJECT_SOURCE_DIR}/source/source_hsolver/module_linalg)
addtest(
  TARGET MODULE_PW_cn_solver
  LIBS parameter base device container psi planewave
  SOURCES test_cn_solver_pw.cpp
          ../cn_solver_pw.cpp
          ../td_orbital_ops.cpp
          ${LINALG_DIR}/hsolver_linear.cpp
          ${LINALG_DIR}/linear_bicgstab.cpp
          ${LINALG_DIR}/linear_cgs.cpp
          ${LINALG_DIR}/linear_gmres.cpp
          ${LINALG_DIR}/linear_algebra.cpp
          ${LINALG_DIR}/cn_subspace.cpp
          ${LINALG_DIR}/linear_low_rank.cpp
          ${LINALG_DIR}/linear_workspace.cpp
          ${LINALG_DIR}/orb_buffer.cpp
          ${LINALG_DIR}/orthonormal.cpp
          ${LINALG_DIR}/orthonormal_small.cpp)

if(ENABLE_MPI)
  addtest(
    TARGET MODULE_PW_td_observables
    LIBS parameter
         base
         device
         container
         psi
         planewave
         nao_smoothgen
         symmetry
    SOURCES td_test_utils.cpp
            test_td_spin_pw.cpp
            test_td_vel_pw.cpp
            ../td_spin_pw.cpp
            ../td_vel_pw.cpp
            ../td_pw.cpp
            ../../module_pwdft/op_pw_vel.cpp
            ../../module_pwdft/projector_gradient.cpp
            ../../module_pwdft/velocity_workspace.cpp
            ../../module_pwdft/nonlocal_workspace.cpp
            ../../module_pwdft/vnl_pw.cpp
            ../../module_pwdft/vnl_pw_init_vnl.cpp
            ../../module_pwdft/vnl_pw_getvnl.cpp
            ../../module_pwdft/vnl_pw_deeq.cpp
            ../../module_pwdft/vnl_pw_grad.cpp
            ../../module_pwdft/vnl_pw_table.cpp
            ../../module_pwdft/vnl_pw_qrad.cpp
            ../../../source_cell/stru_fac.cpp
            ../../../source_cell/stru_fac_k.cpp
            ../../../source_cell/soc.cpp
            ../../module_pwdft/nonlocal_maths.cpp
            ../../module_proj/onsite_proj.cpp
            ../../module_proj/onsite_proj_init.cpp
            ../../module_proj/onsite_proj_overlap.cpp
            ../../module_proj/onsite_proj_tools.cpp
            ../../module_proj/onsite_proj_tools_becp.cpp
            ../../module_proj/radial_proj.cpp
            ../../module_proj/onsite_proj_print.cpp
            ../../../source_basis/module_ao/orb_gaunt_table.cpp
            ../../../source_cell/klist.cpp
            ../../../source_cell/reciprocal_grid.cpp
            ../../../source_cell/klist_io.cpp
            ../../../source_cell/unitcell.cpp
            ../../../source_cell/unitcell_index.cpp
            ../../../source_cell/atom_spec.cpp
            ../../../source_cell/atom_pseudo.cpp
            ../../../source_cell/pseudo.cpp
            ../../../source_cell/magnetism.cpp
            ../../../source_cell/sep.cpp
            ../../../source_cell/sep_cell.cpp
            ../../../source_cell/read_pp.cpp
            ../../../source_cell/read_pp_complete.cpp
            ../../../source_cell/read_pp_upf100.cpp
            ../../../source_cell/read_pp_upf201.cpp
            ../../../source_cell/read_pp_vwr.cpp
            ../../../source_cell/read_pp_blps.cpp
            ../../../source_estate/elecstate.cpp
            ../../../source_estate/module_charge/charge.cpp
            ../../../source_estate/module_pot/potential_new.cpp
            ../../../source_hamilt/module_xc/xc_functional.cpp)
  target_compile_definitions(
    MODULE_PW_td_observables
    PRIVATE TD_TEST_PP_DIR="${PROJECT_SOURCE_DIR}/tests/PP_ORB")
  # Keep real PW initialization and observable implementations without linking the SCF driver.
  if(MSVC)
    target_compile_options(MODULE_PW_td_observables PRIVATE /Gy /Gw)
    target_link_options(MODULE_PW_td_observables PRIVATE /OPT:REF)
  else()
    target_compile_options(MODULE_PW_td_observables PRIVATE -ffunction-sections
                                                            -fdata-sections)
    if(APPLE)
      target_link_options(MODULE_PW_td_observables PRIVATE -Wl,-dead_strip)
    else()
      target_link_options(MODULE_PW_td_observables PRIVATE -Wl,--gc-sections)
    endif()
  endif()
endif()
