add_library(sccs OBJECT sccs_cavity.cpp sccs_parameters.cpp)

if(ENABLE_COVERAGE)
  add_coverage(sccs)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
