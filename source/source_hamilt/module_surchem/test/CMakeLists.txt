abacus_disable_feature_definitions(__LCAO)
abacus_disable_feature_definitions(__LIBXC)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND depend_files 

)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_epsilon
  LIBS parameter planewave device base
  SOURCES cal_epsilon_test.cpp  ../cal_epsilon.cpp ../surchem.cpp 
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_pseudo  
  LIBS parameter planewave device base psi
  SOURCES cal_pseudo_test.cpp  ../cal_pseudo.cpp ../surchem.cpp ../cal_epsilon.cpp ../../../source_pw/module_pwdft/stru_fac.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_totn
  LIBS parameter planewave device base
  SOURCES cal_totn_test.cpp  ../cal_totn.cpp ../surchem.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_vcav
  LIBS parameter planewave device base container
  SOURCES cal_vcav_test.cpp  ../cal_vcav.cpp ../surchem.cpp
  ../../module_xc/xc_grad.cpp ../../module_xc/xc_grad_prepare.cpp ../../module_xc/xc_grad_kernel.cpp ../../module_xc/xc_grad_assemble.cpp ../../module_xc/xc_grad_wfc.cpp ../../module_xc/xc_grad_utils.cpp ../../module_xc/xc_functional.cpp
  ../../module_xc/xc_lda_wrap.cpp ../../module_xc/xc_gga_wrap.cpp
  ../../module_xc/libxc_setup.cpp
  ../../module_xc/libxc_pot.cpp
  ../../module_xc/libxc_lda_wrap.cpp
  ../../module_xc/libxc_gga_wrap.cpp
  ../../module_xc/libxc_mgga_wrap.cpp
  ../../module_xc/xc_gga_corr.cpp ../../module_xc/xc_lda_corr.cpp ../../module_xc/xc_gga_exch.cpp
  ../../module_xc/xc_lda_exch.cpp ../../module_xc/xc_hcth.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_vel
  LIBS parameter planewave device base container
  SOURCES cal_vel_test.cpp  ../cal_vel.cpp ../surchem.cpp ../cal_epsilon.cpp ../minimize_cg.cpp
  ../../module_xc/xc_grad.cpp ../../module_xc/xc_grad_prepare.cpp ../../module_xc/xc_grad_kernel.cpp ../../module_xc/xc_grad_assemble.cpp ../../module_xc/xc_grad_wfc.cpp ../../module_xc/xc_grad_utils.cpp ../../module_xc/xc_functional.cpp
  ../../module_xc/xc_lda_wrap.cpp ../../module_xc/xc_gga_wrap.cpp
  ../../module_xc/libxc_setup.cpp
  ../../module_xc/libxc_pot.cpp
  ../../module_xc/libxc_lda_wrap.cpp
  ../../module_xc/libxc_gga_wrap.cpp
  ../../module_xc/libxc_mgga_wrap.cpp  
  ../../module_xc/xc_gga_corr.cpp ../../module_xc/xc_lda_corr.cpp ../../module_xc/xc_gga_exch.cpp
  ../../module_xc/xc_lda_exch.cpp ../../module_xc/xc_hcth.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_cavity
  SOURCES test_sccs_cavity.cpp ../sccs/sccs_cavity.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pcc
  LIBS base device
  SOURCES test_sccs_pcc.cpp ../pcc/sccs_pcc.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pcc_2d
  LIBS base device
  SOURCES test_sccs_pcc_2d.cpp ../pcc/sccs_pcc.cpp ../pcc/sccs_pcc_2d.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_nonel
  SOURCES test_sccs_nonel.cpp ../sccs/sccs_cavity.cpp ../sccs/sccs_nonel.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_parameters
  SOURCES test_sccs_parameters.cpp ../sccs/sccs_cavity.cpp ../sccs/sccs_parameters.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_poisson
  SOURCES test_sccs_poisson.cpp ../sccs/sccs_poisson.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pw_coulomb
  LIBS planewave device base
  SOURCES test_sccs_pw_coulomb.cpp ../sccs/sccs_pw_coulomb.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_periodic
  LIBS planewave device base
  SOURCES test_sccs_periodic.cpp ../sccs/sccs_cavity.cpp ../sccs/sccs_poisson.cpp
          ../sccs/sccs_pw_coulomb.cpp ../sccs/sccs_pw_reduction.cpp ../sccs/sccs_periodic.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_charge
  SOURCES test_sccs_charge.cpp ../sccs/sccs_charge.cpp ../pcc/sccs_pcc.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pw_charge
  LIBS planewave device base
  SOURCES test_sccs_pw_charge.cpp ../sccs/sccs_pw_charge.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pcc_coulomb
  LIBS planewave device base
  SOURCES test_sccs_pcc_coulomb.cpp ../sccs/sccs_charge.cpp ../pcc/sccs_pcc.cpp
          ../sccs/sccs_poisson.cpp ../sccs/sccs_pw_coulomb.cpp ../pcc/sccs_pcc_coulomb.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pcc_2d_coulomb
  LIBS planewave device base
  SOURCES test_sccs_pcc_2d_coulomb.cpp ../sccs/sccs_cavity.cpp ../sccs/sccs_charge.cpp
          ../sccs/sccs_functional.cpp ../pcc/sccs_pcc.cpp ../pcc/sccs_pcc_2d.cpp
          ../sccs/sccs_poisson.cpp ../sccs/sccs_periodic.cpp ../sccs/sccs_pw_coulomb.cpp
          ../pcc/sccs_pcc_2d_coulomb.cpp ../sccs/sccs_pw_charge.cpp
          ../sccs/sccs_pw_reduction.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_functional
  LIBS planewave device base
  SOURCES test_sccs_functional.cpp ../sccs/sccs_charge.cpp ../pcc/sccs_pcc.cpp
          ../sccs/sccs_functional.cpp ../sccs/sccs_pw_coulomb.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pw_force
  LIBS planewave device base cell symmetry neighbor_search neighbor
  SOURCES test_sccs_pw_force.cpp ../sccs/sccs_pw_charge.cpp ../sccs/sccs_pw_force.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_pw_nonel
  LIBS planewave device base
  SOURCES test_sccs_pw_nonel.cpp ../sccs/sccs_charge.cpp ../pcc/sccs_pcc.cpp
          ../sccs/sccs_pw_coulomb.cpp ../sccs/sccs_pw_nonel.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sccs_driver
  LIBS planewave device base
  SOURCES test_sccs_driver.cpp ../sccs/sccs_cavity.cpp ../sccs/sccs_charge.cpp
          ../sccs/sccs_driver.cpp ../sccs/sccs_functional.cpp ../sccs/sccs_nonel.cpp
          ../sccs/sccs_parameters.cpp ../pcc/sccs_pcc.cpp ../pcc/sccs_pcc_2d.cpp
          ../pcc/sccs_pcc_2d_coulomb.cpp ../pcc/sccs_pcc_coulomb.cpp
          ../sccs/sccs_periodic.cpp ../sccs/sccs_poisson.cpp ../sccs/sccs_pw_charge.cpp
          ../sccs/sccs_pw_coulomb.cpp ../sccs/sccs_pw_nonel.cpp ../sccs/sccs_pw_reduction.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_h_corr_sccs
  LIBS planewave device base cell symmetry neighbor_search neighbor
  SOURCES test_h_corr_sccs.cpp test_h_corr_pcc.cpp test_surchem_input.cpp ../surchem_input.cpp ../surchem.cpp ../sccs/h_corr_sccs.cpp ../pcc/h_corr_pcc.cpp
          ../sccs/sccs_cavity.cpp ../sccs/sccs_charge.cpp ../sccs/sccs_driver.cpp
          ../sccs/sccs_functional.cpp ../sccs/sccs_nonel.cpp ../sccs/sccs_parameters.cpp
          ../pcc/sccs_pcc.cpp ../pcc/sccs_pcc_2d.cpp ../pcc/sccs_pcc_2d_coulomb.cpp
          ../pcc/sccs_pcc_coulomb.cpp ../sccs/sccs_periodic.cpp
          ../sccs/sccs_poisson.cpp ../sccs/sccs_pw_charge.cpp ../sccs/sccs_pw_coulomb.cpp
          ../sccs/sccs_pw_force.cpp ../sccs/sccs_pw_nonel.cpp ../sccs/sccs_pw_reduction.cpp
)

AddTest(
  TARGET MODULE_HAMILT_surchem_sol_force
  LIBS planewave device base cell symmetry neighbor_search neighbor
  SOURCES test_sol_force.cpp ../surchem.cpp ../sccs/h_corr_sccs.cpp ../pcc/h_corr_pcc.cpp ../sol_force.cpp
          ../sccs/sccs_cavity.cpp ../sccs/sccs_charge.cpp ../sccs/sccs_driver.cpp
          ../sccs/sccs_functional.cpp ../sccs/sccs_nonel.cpp ../sccs/sccs_parameters.cpp
          ../pcc/sccs_pcc.cpp ../pcc/sccs_pcc_2d.cpp ../pcc/sccs_pcc_2d_coulomb.cpp
          ../pcc/sccs_pcc_coulomb.cpp ../sccs/sccs_periodic.cpp
          ../sccs/sccs_poisson.cpp ../sccs/sccs_pw_charge.cpp ../sccs/sccs_pw_coulomb.cpp
          ../sccs/sccs_pw_force.cpp ../sccs/sccs_pw_nonel.cpp ../sccs/sccs_pw_reduction.cpp
)
