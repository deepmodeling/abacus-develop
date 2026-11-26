abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

addtest(
  TARGET MODULE_HAMILT_vdwd3EvaluatorTest
  LIBS parameter base device vdw
  SOURCES vdwd3_evaluator_test.cpp)
