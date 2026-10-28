abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

addtest(
  TARGET
  MODULE_RPA_librpa_stru_symmetry_test
  LIBS
  base
  device
  symmetry
  SOURCES
  test_librpa_stru_symmetry.cpp
  ../librpa_bz_sampling.cpp
  ../librpa_stru_symmetry.cpp
  ../librpa_stru_units.cpp)

add_executable(MODULE_RPA_librpa_stru_units_test test_librpa_stru_units.cpp
                                                 ../librpa_stru_units.cpp)
target_include_directories(MODULE_RPA_librpa_stru_units_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
install(TARGETS MODULE_RPA_librpa_stru_units_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_stru_units_test
  COMMAND ./MODULE_RPA_librpa_stru_units_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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

add_executable(MODULE_RPA_librpa_stru_writers_test test_librpa_stru_writers.cpp
                                                   ../librpa_stru_symmetry.cpp)
target_include_directories(MODULE_RPA_librpa_stru_writers_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RPA_librpa_stru_writers_test
                      PRIVATE base device symmetry abacus::linalg_libs)
install(TARGETS MODULE_RPA_librpa_stru_writers_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_stru_writers_test
  COMMAND ./MODULE_RPA_librpa_stru_writers_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
