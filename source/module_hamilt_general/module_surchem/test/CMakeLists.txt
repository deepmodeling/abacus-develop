remove_definitions(-D__LCAO )
remove_definitions(-DUSE_LIBXC)
remove_definitions(-DUSE_PAW)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND depend_files 

)

AddTest(
  TARGET surchem_cal_epsilon
  LIBS ${math_libs} planewave device base
  SOURCES cal_epsilon_test.cpp  ../cal_epsilon.cpp ../surchem.cpp 
)

AddTest(
  TARGET surchem_cal_pseudo  LIBS ${math_libs} planewave device base psi
  SOURCES cal_pseudo_test.cpp  ../cal_pseudo.cpp ../surchem.cpp ../cal_epsilon.cpp ../../../module_hamilt_pw/hamilt_pwdft/structure_factor.cpp
  ../../../module_hamilt_pw/hamilt_pwdft/parallel_grid.cpp 

)

AddTest(
  TARGET surchem_cal_totn
  LIBS ${math_libs} planewave device base
  SOURCES cal_totn_test.cpp  ../cal_totn.cpp ../surchem.cpp ../../../module_hamilt_pw/hamilt_pwdft/parallel_grid.cpp 
)

AddTest(
  TARGET surchem_cal_vcav
  LIBS ${math_libs} planewave device base
  SOURCES cal_vcav_test.cpp  ../cal_vcav.cpp ../surchem.cpp ../../../module_hamilt_pw/hamilt_pwdft/parallel_grid.cpp 
  ../../module_xc/xc_functional_gradcorr.cpp ../../module_xc/xc_functional.cpp
  ../../module_xc/xc_functional_wrapper_xc.cpp ../../module_xc/xc_functional_wrapper_gcxc.cpp
  ../../module_xc/xc_functional_wrapper_tauxc.cpp
  ../../module_xc/xc_funct_corr_gga.cpp ../../module_xc/xc_funct_corr_lda.cpp ../../module_xc/xc_funct_exch_gga.cpp
  ../../module_xc/xc_funct_exch_lda.cpp ../../module_xc/xc_funct_hcth.cpp
)

AddTest(
  TARGET surchem_cal_vel
  LIBS ${math_libs} planewave device base
  SOURCES cal_vel_test.cpp  ../cal_vel.cpp ../surchem.cpp ../cal_epsilon.cpp ../minimize_cg.cpp ../../../module_hamilt_pw/hamilt_pwdft/parallel_grid.cpp 
  ../../module_xc/xc_functional_gradcorr.cpp ../../module_xc/xc_functional.cpp
  ../../module_xc/xc_functional_wrapper_xc.cpp ../../module_xc/xc_functional_wrapper_gcxc.cpp
  ../../module_xc/xc_functional_wrapper_tauxc.cpp
  ../../module_xc/xc_funct_corr_gga.cpp ../../module_xc/xc_funct_corr_lda.cpp ../../module_xc/xc_funct_exch_gga.cpp
  ../../module_xc/xc_funct_exch_lda.cpp ../../module_xc/xc_funct_hcth.cpp
)