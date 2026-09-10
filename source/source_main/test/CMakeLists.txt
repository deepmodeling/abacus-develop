abacus_disable_feature_definitions(__MPI)

AddTest(
  TARGET MODULE_MAIN_para_worlds_global
  SOURCES test_para_worlds_global.cpp ../para_worlds_global.cpp
  LIBS base device
)

# The sources use repository-root include paths (e.g. "source_base/...");
# OBJECT library `base` does not propagate its include interface here, so add
# the source root explicitly.
target_include_directories(MODULE_MAIN_para_worlds_global PRIVATE
  ${CMAKE_SOURCE_DIR}/source)
