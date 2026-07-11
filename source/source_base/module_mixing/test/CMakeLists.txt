abacus_disable_feature_definitions(__MPI)
AddTest(
  TARGET MODULE_BASE_MIXING_unittests
  LIBS parameter base device
  SOURCES mixing_test.cpp
)