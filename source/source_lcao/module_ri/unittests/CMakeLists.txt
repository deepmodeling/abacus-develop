abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
if(ENABLE_LCAO)
  AddTest(
    TARGET MODULE_RI_exx_workspace
    LIBS base parameter device ${math_libs}
    SOURCES test_exx_workspace.cpp ../exx_workspace.cpp ../../../source_basis/module_ao/orb_atomic_lm.cpp
  )
endif()
AddTest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES test_mix_dmk_2d.cpp ../mix_dmk_2d.cpp
)
AddTest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES test_write_ri_cv.cpp
)
AddTest(
  TARGET MODULE_RI_abfs_vector3_order
  LIBS parameter
  SOURCES test_abfs_vector3_order.cpp
)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/support"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
add_executable(MODULE_RI_ewald_tail_utils_test test_ewald_tail_utils.cpp)
install(TARGETS MODULE_RI_ewald_tail_utils_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_tail_utils_test
  COMMAND ./MODULE_RI_ewald_tail_utils_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_cutoff_test test_gaussian_abfs_cutoff.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_cutoff_test
  COMMAND ./MODULE_RI_gaussian_abfs_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_gaussian_abfs_2d_integrals_test test_gaussian_abfs_2d_integrals.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_2d_integrals_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_2d_integrals_test
  COMMAND ./MODULE_RI_gaussian_abfs_2d_integrals_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(MODULE_RI_ewald_gaussian_cutoff_test test_ewald_gaussian_cutoff.cpp)
install(TARGETS MODULE_RI_ewald_gaussian_cutoff_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_gaussian_cutoff_test
  COMMAND ./MODULE_RI_ewald_gaussian_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_2d_coulomb_test
  test_direct_2d_coulomb.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_2d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_2d_coulomb_test PRIVATE base device parameter abacus::linalg_libs)
install(TARGETS MODULE_RI_direct_2d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_2d_coulomb_test
  COMMAND ./MODULE_RI_direct_2d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_executable(
  MODULE_RI_direct_3d_coulomb_test
  test_direct_3d_coulomb.cpp
  ../direct_2d_coulomb.cpp
)
target_include_directories(MODULE_RI_direct_3d_coulomb_test PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_3d_coulomb_test PRIVATE base device parameter abacus::linalg_libs)
install(TARGETS MODULE_RI_direct_3d_coulomb_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_3d_coulomb_test
  COMMAND ./MODULE_RI_direct_3d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
if(ENABLE_MPI)
  add_executable(MODULE_RI_ewald_mpi_distribution_test test_ewald_mpi_distribution.cpp)
  target_link_libraries(MODULE_RI_ewald_mpi_distribution_test PRIVATE MPI::MPI_CXX)
  install(TARGETS MODULE_RI_ewald_mpi_distribution_test DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(
    NAME MODULE_RI_ewald_mpi_distribution_parallel_test
    COMMAND mpirun -np 4 ./MODULE_RI_ewald_mpi_distribution_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
