remove_definitions(-D__MPI)
list(APPEND depend_files 
  ../../../module_base/global_variable.cpp
  
)
AddTest(
  TARGET surchem_atom_in
  SOURCES atom_in_test.cpp
)
AddTest(
  TARGET surchem_cal_epsilon
  SOURCES cal_epsilon_test.cpp  ../cal_epsilon.cpp
  ${depend_files}
)
