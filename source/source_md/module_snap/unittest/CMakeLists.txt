addtest(TARGET MODULE_MD_SNAP_model SOURCES test_snap_model.cpp
                                            ../snap_model.cpp)

addtest(TARGET MODULE_MD_SNAP_cg SOURCES test_snap_cg.cpp ../snap_cg.cpp)

addtest(
  TARGET MODULE_MD_SNAP_descriptor
  LIBS base device
  SOURCES test_snap_descriptor.cpp ../snap_descriptor.cpp
          ../snap_derivative.cpp ../snap_model.cpp ../snap_cg.cpp)

addtest(
  TARGET MODULE_MD_SNAP_calculator
  LIBS base device
  SOURCES test_snap_calculator.cpp ../snap_calculator.cpp
          ../snap_descriptor.cpp ../snap_derivative.cpp ../snap_model.cpp
          ../snap_cg.cpp)
