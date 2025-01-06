remove_definitions(-DUSE_LIBXC)
AddTest(
  TARGET AX_test
  LIBS parameter base ${math_libs} container device psi 
  SOURCES AX_test.cpp ../../utils/lr_util.cpp ../AX_parallel.cpp ../AX_serial.cpp
)