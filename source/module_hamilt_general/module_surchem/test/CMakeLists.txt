remove_definitions(-D__LCAO )
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND depend_files 

)
AddTest(
  TARGET surchem_atom_in
  SOURCES atom_in_test.cpp
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