addtest(
  TARGET MODULE_RPA_rpa_lri_detail_test
  LIBS parameter base device
  SOURCES test_rpa_lri_detail.cpp ../rpa_lri_detail.cpp ../rpa_lri_coeff_v1.cpp
          ../rpa_lri_bz_io.cpp ../module_librpa_io/librpa_bz_sampling.cpp)
addtest(
  TARGET MODULE_RPA_rpa_lri_band_io_test
  LIBS parameter base device
  SOURCES test_rpa_lri_band_io.cpp ../rpa_lri_band_io.cpp)
