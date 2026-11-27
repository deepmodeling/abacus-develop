addtest(
  TARGET MODULE_PW_radial_proj_test
  LIBS parameter base device
  SOURCES test_radial_proj.cpp ../radial_proj.cpp)

addtest(
  TARGET MODULE_PW_projgen_test
  LIBS parameter base device
  SOURCES test_projgen.cpp ../projgen.cpp)
