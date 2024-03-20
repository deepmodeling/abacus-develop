remove_definitions(-D __MPI)

add_library(tddft_test_lib tddft_test.cpp)

AddTest(
  TARGET tddft_middle_hamilt_test
  LIBS ${math_libs} base device tddft_test_lib  
  SOURCES middle_hamilt_test.cpp ../middle_hamilt.cpp 
)

AddTest(
  TARGET tddft_bandenergy_test
  LIBS ${math_libs} base device tddft_test_lib  
  SOURCES bandenergy_test.cpp ../bandenergy.cpp 
)

AddTest(
  TARGET tddft_norm_psi_test
  LIBS ${math_libs} base device tddft_test_lib  
  SOURCES norm_psi_test.cpp ../norm_psi.cpp 
)

AddTest(
  TARGET tddft_upsi_test
  LIBS ${math_libs} base device tddft_test_lib  
  SOURCES upsi_test1.cpp upsi_test2.cpp upsi_test3.cpp ../upsi.cpp 
)

AddTest(
  TARGET tddft_propagator_test
  LIBS ${math_libs} base device tddft_test_lib  
  SOURCES propagator_test1.cpp propagator_test2.cpp propagator_test3.cpp ../propagator.cpp 
)

