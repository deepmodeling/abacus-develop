remove_definitions(-D__MPI)
AddTest(
  TARGET test_mixing
  LIBS parameter base device ${math_libs} 
  SOURCES mixing_test.cpp
)