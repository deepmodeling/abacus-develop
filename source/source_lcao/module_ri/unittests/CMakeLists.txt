abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
if(ENABLE_LCAO)
  addtest(
    TARGET
    MODULE_RI_exx_workspace
    LIBS
    base
    parameter
    device
    ${math_libs}
    SOURCES
    test_exx_workspace.cpp
    ../exx_workspace.cpp
    ../lri_cv.cpp
    ../../../source_basis/module_ao/orb_atomic_lm.cpp)
  # Link the compiled LRI lifecycle while keeping unrelated integral assembly out
  # of this focused geometry-sharing test.
  target_compile_options(MODULE_RI_exx_workspace PRIVATE -ffunction-sections -fdata-sections)
  target_link_options(MODULE_RI_exx_workspace PRIVATE -Wl,--gc-sections)
endif()
addtest(
  TARGET
  MODULE_RI_dm_mixing_test
  LIBS
  parameter
  base
  device
  SOURCES
  test_mix_dmk_2d.cpp
  ../mix_dmk_2d.cpp)
addtest(
  TARGET
  MODULE_RI_ri_cv_io_test
  LIBS
  base
  device
  parameter
  SOURCES
  test_write_ri_cv.cpp)
addtest(TARGET MODULE_RI_abfs_vector3_order LIBS parameter SOURCES
        test_abfs_vector3_order.cpp)
abacus_copy_test_files(support)

add_executable(MODULE_RI_ewald_tail_utils_test test_ewald_tail_utils.cpp)
install(TARGETS MODULE_RI_ewald_tail_utils_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_tail_utils_test
  COMMAND ./MODULE_RI_ewald_tail_utils_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RI_gaussian_abfs_cutoff_test
               test_gaussian_abfs_cutoff.cpp)
install(TARGETS MODULE_RI_gaussian_abfs_cutoff_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_abfs_cutoff_test
  COMMAND ./MODULE_RI_gaussian_abfs_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RI_gaussian_2d_integrals_test
               test_gaussian_2d_integrals.cpp)
install(TARGETS MODULE_RI_gaussian_2d_integrals_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_2d_integrals_test
  COMMAND ./MODULE_RI_gaussian_2d_integrals_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RI_ewald_cutoff_test test_ewald_cutoff.cpp)
install(TARGETS MODULE_RI_ewald_cutoff_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_ewald_cutoff_test
  COMMAND ./MODULE_RI_ewald_cutoff_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(
  MODULE_RI_gaussian_dipole_test
  test_gaussian_dipole.cpp ../gaussian_abfs.cpp
  ../../../source_basis/module_ao/orb_atomic_lm.cpp)
target_include_directories(MODULE_RI_gaussian_dipole_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(
  MODULE_RI_gaussian_dipole_test
  PRIVATE base
          device
          parameter
          cell
          symmetry
          neighbor_search
          domain_decomposition
          neighbor
          abacus::linalg_libs)
install(TARGETS MODULE_RI_gaussian_dipole_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_gaussian_dipole_test
  COMMAND ./MODULE_RI_gaussian_dipole_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RI_direct_2d_coulomb_test test_direct_2d_coulomb.cpp
                                                ../direct_2d_coulomb.cpp)
target_include_directories(MODULE_RI_direct_2d_coulomb_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_2d_coulomb_test
                      PRIVATE base device parameter abacus::linalg_libs)
install(TARGETS MODULE_RI_direct_2d_coulomb_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_2d_coulomb_test
  COMMAND ./MODULE_RI_direct_2d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_executable(MODULE_RI_direct_3d_coulomb_test test_direct_3d_coulomb.cpp
                                                ../direct_2d_coulomb.cpp)
target_include_directories(MODULE_RI_direct_3d_coulomb_test
                           PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(MODULE_RI_direct_3d_coulomb_test
                      PRIVATE base device parameter abacus::linalg_libs)
install(TARGETS MODULE_RI_direct_3d_coulomb_test
        DESTINATION ${CMAKE_BINARY_DIR}/tests)
add_test(
  NAME MODULE_RI_direct_3d_coulomb_test
  COMMAND ./MODULE_RI_direct_3d_coulomb_test
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(ENABLE_MPI)
  add_executable(MODULE_RI_ewald_mpi_dist_test test_ewald_mpi_dist.cpp)
  target_link_libraries(
    MODULE_RI_ewald_mpi_dist_test PRIVATE MPI::MPI_CXX base device parameter
                                          abacus::linalg_libs)
  install(TARGETS MODULE_RI_ewald_mpi_dist_test
          DESTINATION ${CMAKE_BINARY_DIR}/tests)
  add_test(
    NAME MODULE_RI_ewald_mpi_dist_parallel_test
    COMMAND mpirun -np 4 "$<TARGET_FILE:MODULE_RI_ewald_mpi_dist_test>"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

addtest(
  TARGET MODULE_RI_invmat_ri
  LIBS parameter base device ${math_libs}
  SOURCES test_invmat_ri.cpp ../invmat_ri.cpp ../lri_cv_tools.cpp
          ../abfs_construct_pca.cpp
          ../../../source_basis/module_ao/orb_atomic_lm.cpp)
# Keep only the tensor eigensolver from the PCA translation unit in this focused test.
target_compile_options(MODULE_RI_invmat_ri PRIVATE -ffunction-sections -fdata-sections)
target_link_options(MODULE_RI_invmat_ri PRIVATE -Wl,--gc-sections)
