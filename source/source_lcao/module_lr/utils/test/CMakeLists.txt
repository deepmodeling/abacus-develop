abacus_disable_feature_definitions(__LIBXC)
AddTest(
  TARGET MODULE_LR_lr_util_phys_test
  LIBS parameter base device container planewave #for FFT
  SOURCES lr_util_physics_test.cpp  ../lr_util.cpp
  
)

AddTest(
  TARGET MODULE_LR_lr_util_algo_test
  LIBS parameter base device psi container planewave #for FFT
  SOURCES lr_util_algorithms_test.cpp ../lr_util.cpp 
)