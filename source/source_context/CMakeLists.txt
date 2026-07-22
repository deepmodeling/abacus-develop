add_library(
  context
  OBJECT
  simulation_context_builder.cpp
  simulation_context_binding.cpp)

target_compile_definitions(
  context
  PRIVATE
  ABACUS_CAN_BIND_SIMULATION_CONTEXT
  ABACUS_CAN_READ_SIMULATION_CONTEXT)

if(ENABLE_COVERAGE)
  add_coverage(context)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
