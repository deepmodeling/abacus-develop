AddTest(
  TARGET MODULE_MD_SNAP_model
  SOURCES test_model.cpp ../model.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_clebsch_gordan
  SOURCES test_cg.cpp ../cg.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_descriptor
  LIBS base device
  SOURCES test_desc.cpp ../desc.cpp ../deriv.cpp ../model.cpp ../cg.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_calculator
  LIBS base device
  SOURCES test_calc.cpp ../calc.cpp ../desc.cpp ../deriv.cpp ../model.cpp ../cg.cpp
)
