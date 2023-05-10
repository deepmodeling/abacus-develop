install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND depend_files 
  ../../../module_base/global_variable.cpp
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