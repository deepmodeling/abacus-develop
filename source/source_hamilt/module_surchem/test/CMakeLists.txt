remove_definitions(-D__LCAO )
remove_definitions(-DUSE_LIBXC)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND depend_files 

)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_epsilon
  LIBS parameter ${math_libs} planewave device base 
  SOURCES cal_epsilon_test.cpp  ../cal_epsilon.cpp ../surchem.cpp 
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_pseudo  
  LIBS parameter ${math_libs} planewave device base psi 
  SOURCES cal_pseudo_test.cpp  ../cal_pseudo.cpp ../surchem.cpp ../cal_epsilon.cpp ../../../source_pw/module_pwdft/structure_factor.cpp
  ../../../source_pw/module_pwdft/parallel_grid.cpp 

)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_totn
  LIBS parameter ${math_libs} planewave device base 
  SOURCES cal_totn_test.cpp  ../cal_totn.cpp ../surchem.cpp ../../../source_pw/module_pwdft/parallel_grid.cpp 
)

AddTest(
  TARGET MODULE_HAMILT_surchem_cal_vcav
  LIBS parameter ${math_libs} planewave device base container 
  SOURCES cal_vcav_test.cpp  ../cal_vcav.cpp ../surchem.cpp ../../../source_pw/module_pwdft/parallel_grid.cpp 
  ../../module_xc/xc_grad.cpp ../../module_xc/xc_functional.cpp
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
  LIBS parameter ${math_libs} planewave device base container 
  SOURCES cal_vel_test.cpp  ../cal_vel.cpp ../surchem.cpp ../cal_epsilon.cpp ../minimize_cg.cpp ../../../source_pw/module_pwdft/parallel_grid.cpp 
  ../../module_xc/xc_grad.cpp ../../module_xc/xc_functional.cpp
  ../../module_xc/xc_lda_wrap.cpp ../../module_xc/xc_gga_wrap.cpp
  ../../module_xc/libxc_setup.cpp
  ../../module_xc/libxc_pot.cpp
  ../../module_xc/libxc_lda_wrap.cpp
  ../../module_xc/libxc_gga_wrap.cpp
  ../../module_xc/libxc_mgga_wrap.cpp  
  ../../module_xc/xc_gga_corr.cpp ../../module_xc/xc_lda_corr.cpp ../../module_xc/xc_gga_exch.cpp
  ../../module_xc/xc_lda_exch.cpp ../../module_xc/xc_hcth.cpp
)