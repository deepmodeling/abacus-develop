abacus_disable_feature_definitions(__MPI)
AddTest(
  TARGET MODULE_BASE_MIXING_unittests
  LIBS parameter base device ${math_libs} 
  SOURCES mixing_test.cpp
)