abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES dm_mixing_test.cpp ../mix_dmk_2d.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES ri_cv_io_test.cpp
)
add_executable(MODULE_RI_ewald_tail_utils_test ewald_tail_utils_test.cpp)
install(TARGETS MODULE_RI_ewald_tail_utils_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_tail_utils_test
  COMMAND ./MODULE_RI_ewald_tail_utils_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_cutoff_test gaussian_abfs_cutoff_test.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_cutoff_test
  COMMAND ./MODULE_RI_gaussian_abfs_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_2d_integrals_test gaussian_abfs_2d_integrals_test.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_2d_integrals_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_2d_integrals_test
  COMMAND ./MODULE_RI_gaussian_abfs_2d_integrals_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_ewald_gaussian_cutoff_test ewald_gaussian_cutoff_test.cpp)
install(TARGETS MODULE_RI_ewald_gaussian_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_gaussian_cutoff_test
  COMMAND ./MODULE_RI_ewald_gaussian_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_2d_coulomb_test
  direct_2d_coulomb_test.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_2d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_2d_coulomb_test PRIVATE base ${math_libs} device parameter)
install(TARGETS MODULE_RI_direct_2d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_2d_coulomb_test
  COMMAND ./MODULE_RI_direct_2d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_3d_coulomb_test
  direct_3d_coulomb_test.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_3d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_3d_coulomb_test PRIVATE base ${math_libs} device parameter)
install(TARGETS MODULE_RI_direct_3d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_3d_coulomb_test
  COMMAND ./MODULE_RI_direct_3d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
if(ENABLE_MPI)
  add_executable(MODULE_RI_ewald_mpi_distribution_test ewald_mpi_distribution_test.cpp)
  target_link_libraries(MODULE_RI_ewald_mpi_distribution_test PRIVATE MPI::MPI_CXX)
  install(TARGETS MODULE_RI_ewald_mpi_distribution_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(
    NAME MODULE_RI_ewald_mpi_distribution_test_parallel
    COMMAND mpirun -np 4 ./MODULE_RI_ewald_mpi_distribution_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
