abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
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
AddTest(
  TARGET MODULE_RI_rpa_lri_detail_test
  LIBS parameter base device
  SOURCES test_rpa_lri_detail.cpp ../rpa_lri_detail.cpp ../rpa_lri_coeff_v1.cpp ../rpa_lri_bz_io.cpp
)
AddTest(
  TARGET MODULE_RI_rpa_lri_band_io_test
  LIBS parameter base device
  SOURCES test_rpa_lri_band_io.cpp ../rpa_lri_band_io.cpp
)
AddTest(
  TARGET MODULE_RI_rpa_lri_symmetry_test
  LIBS base device symmetry
  SOURCES test_rpa_lri_symmetry.cpp
)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/support"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
