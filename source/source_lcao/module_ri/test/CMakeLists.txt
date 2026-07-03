abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base ${math_libs} device
  SOURCES dm_mixing_test.cpp ../Mix_DMk_2D.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base ${math_libs} device parameter
  SOURCES ri_cv_io_test.cpp
)
AddTest(
  TARGET MODULE_RI_abfs-vector3_order
  LIBS parameter
  SOURCES abfs-vector3_order_test.cpp
)
install(DIRECTORY support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
