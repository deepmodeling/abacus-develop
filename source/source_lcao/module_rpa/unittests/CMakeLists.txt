AddTest(
  TARGET MODULE_RPA_rpa_lri_detail_test
  LIBS parameter base device
  SOURCES test_rpa_lri_detail.cpp
          ../rpa_lri_detail.cpp
          ../rpa_lri_coeff_v1.cpp
          ../rpa_lri_bz_io.cpp
          ../module_librpa_io/librpa_bz_sampling.cpp
)
AddTest(
  TARGET MODULE_RPA_rpa_lri_band_io_test
  LIBS parameter base device
  SOURCES test_rpa_lri_band_io.cpp ../rpa_lri_band_io.cpp
)

AddTest(
  TARGET MODULE_RPA_rpa_abfs_preorthogonalization_test
  LIBS parameter base ${math_libs} device
  SOURCES test_rpa_abfs_preorthogonalization.cpp
          ../rpa_abfs_preorthogonalization.cpp
          ../../../source_basis/module_ao/orb_atomic_lm.cpp
)
