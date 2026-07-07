abacus_disable_feature_definitions(__CUDA)

AddTest(
  TARGET dftu_pw_test
  LIBS base device parameter
  SOURCES dftu_pw_test.cpp
)

AddTest(
  TARGET dftu_core_test
  LIBS base device
  SOURCES dftu_core_test.cpp
)

AddTest(
  TARGET dftu_operator_test
  LIBS base device
  SOURCES dftu_operator_test.cpp
)
