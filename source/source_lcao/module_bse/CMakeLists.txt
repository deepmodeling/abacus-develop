add_library(
  bse
  OBJECT
  hamilt_bse.cpp
  hamilt_bse_solver.cpp
  bse_util.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(bse)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
