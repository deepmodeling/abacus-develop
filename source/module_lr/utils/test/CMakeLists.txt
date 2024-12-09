remove_definitions(-DUSE_LIBXC)
AddTest(
  TARGET lr_util_phys_test
  LIBS parameter base ${math_libs} device container planewave #for FFT
  SOURCES lr_util_physics_test.cpp  ../lr_util.cpp
  ../../../module_io/orb_io.cpp
)

AddTest(
  TARGET lr_util_algo_test
  LIBS parameter base ${math_libs} device psi container planewave #for FFT
  SOURCES lr_util_algorithms_test.cpp ../lr_util.cpp 
)