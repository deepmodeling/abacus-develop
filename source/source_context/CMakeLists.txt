add_library(
  context
  OBJECT
  simulation_context_builder.cpp)

if(ENABLE_COVERAGE)
  add_coverage(context)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
