AddTest(
  TARGET MODULE_MD_SNAP_model
  SOURCES test_snap_model.cpp ../snap_model.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_clebsch_gordan
  SOURCES test_snap_clebsch_gordan.cpp ../snap_clebsch_gordan.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_descriptor
  LIBS base device
  SOURCES test_snap_descriptor.cpp ../snap_descriptor.cpp ../snap_descriptor_derivative.cpp ../snap_model.cpp ../snap_clebsch_gordan.cpp
)

AddTest(
  TARGET MODULE_MD_SNAP_calculator
  LIBS base device
  SOURCES test_snap_calculator.cpp ../snap_calculator.cpp ../snap_descriptor.cpp ../snap_descriptor_derivative.cpp ../snap_model.cpp ../snap_clebsch_gordan.cpp
)
