abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
addtest(TARGET MODULE_RI_smooth_cut_coulomb SOURCES test_smooth_cut_coulomb.cpp
                                                    ../smooth_cut_coulomb.cpp)
if(ENABLE_LCAO)
  addtest(
    TARGET MODULE_RI_exx_workspace
    LIBS base parameter device ${math_libs}
    SOURCES test_exx_workspace.cpp ../exx_workspace.cpp ../lri_cv.cpp
            ../../../source_basis/module_ao/orb_atomic_lm.cpp)
  # Link the compiled LRI lifecycle while keeping unrelated integral assembly out
  # of this focused geometry-sharing test.
  target_compile_options(MODULE_RI_exx_workspace PRIVATE -ffunction-sections
                                                         -fdata-sections)
  target_link_options(MODULE_RI_exx_workspace PRIVATE -Wl,--gc-sections)
endif()
addtest(
  TARGET MODULE_RI_dm_mixing_test
  LIBS parameter base device
  SOURCES test_mix_dmk_2d.cpp ../mix_dmk_2d.cpp)
addtest(
  TARGET MODULE_RI_ri_cv_io_test
  LIBS base device parameter
  SOURCES test_write_ri_cv.cpp)
addtest(
  TARGET MODULE_RI_abfs_vector3_order
  LIBS parameter
  SOURCES test_abfs_vector3_order.cpp)
abacus_copy_test_files(support)

addtest(
  TARGET MODULE_RI_invmat_ri
  LIBS parameter base device ${math_libs}
  SOURCES test_invmat_ri.cpp ../invmat_ri.cpp ../lri_cv_tools.cpp
          ../abfs_construct_pca.cpp
          ../../../source_basis/module_ao/orb_atomic_lm.cpp)
# Keep only the tensor eigensolver from the PCA translation unit in this focused test.
target_compile_options(MODULE_RI_invmat_ri PRIVATE -ffunction-sections
                                                   -fdata-sections)
target_link_options(MODULE_RI_invmat_ri PRIVATE -Wl,--gc-sections)
