abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

add_executable(MODULE_RPA_librpa_bz_sampling_test test_librpa_bz_sampling.cpp
                                                  ../librpa_bz_sampling.cpp)
target_include_directories(MODULE_RPA_librpa_bz_sampling_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
install(TARGETS MODULE_RPA_librpa_bz_sampling_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_bz_sampling_test
  COMMAND ./MODULE_RPA_librpa_bz_sampling_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RPA_librpa_2d_head_test test_librpa_2d_head.cpp)
target_include_directories(MODULE_RPA_librpa_2d_head_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
install(TARGETS MODULE_RPA_librpa_2d_head_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_2d_head_test
  COMMAND ./MODULE_RPA_librpa_2d_head_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
