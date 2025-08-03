remove_definitions(-DUSE_LIBXC)
AddTest(
  TARGET ao_to_mo_test
  LIBS parameter base ${math_libs} container device psi 
  SOURCES ao_to_mo_test.cpp ../../utils/lr_util.cpp ../ao_to_mo_parallel.cpp ../ao_to_mo_serial.cpp
)