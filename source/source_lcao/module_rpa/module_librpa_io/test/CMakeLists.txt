abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

add_executable(MODULE_RPA_librpa_stru_units_test librpa_stru_units_test.cpp)
install(TARGETS MODULE_RPA_librpa_stru_units_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_stru_units_test
  COMMAND ./MODULE_RPA_librpa_stru_units_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(
  MODULE_RPA_librpa_stru_symmetry_test
  librpa_stru_symmetry_test.cpp
  ${ABACUS_SOURCE_DIR}/source_cell/module_symmetry/symm_rot_spin.cpp
)
target_link_libraries(MODULE_RPA_librpa_stru_symmetry_test
  PRIVATE base device abacus::linalg_libs)
install(TARGETS MODULE_RPA_librpa_stru_symmetry_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_stru_symmetry_test
  COMMAND ./MODULE_RPA_librpa_stru_symmetry_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(MODULE_RPA_librpa_bz_sampling_test librpa_bz_sampling_test.cpp)
install(TARGETS MODULE_RPA_librpa_bz_sampling_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_bz_sampling_test
  COMMAND ./MODULE_RPA_librpa_bz_sampling_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(MODULE_RPA_librpa_2d_coulomb_head_test librpa_2d_coulomb_head_test.cpp)
install(TARGETS MODULE_RPA_librpa_2d_coulomb_head_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RPA_librpa_2d_coulomb_head_test
  COMMAND ./MODULE_RPA_librpa_2d_coulomb_head_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
