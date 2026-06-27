abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)

AddTest(
  TARGET MODULE_IO_JSON_OUTPUT_TEST
  LIBS parameter ${math_libs} base device cell_info json_output
  SOURCES para_json_test.cpp ../para_json.cpp
)
