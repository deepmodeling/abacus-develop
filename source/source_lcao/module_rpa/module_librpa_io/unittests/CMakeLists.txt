addtest(
  TARGET MODULE_RPA_librpa_stru_symmetry_test
  LIBS base device symmetry
  SOURCES test_librpa_stru_symmetry.cpp ../librpa_bz_sampling.cpp
          ../librpa_stru_symmetry.cpp ../librpa_stru_units.cpp)
