AddTest(
  TARGET MODULE_CONTEXT_builder_test
  LIBS context parameter base device
  SOURCES simulation_context_builder_test.cpp
)

target_compile_definitions(
  MODULE_CONTEXT_builder_test
  PRIVATE
  ABACUS_CAN_BIND_SIMULATION_CONTEXT
  ABACUS_CAN_READ_SIMULATION_CONTEXT)

AddTest(
  TARGET MODULE_CONTEXT_binding_test
  LIBS context parameter base device
  SOURCES simulation_context_binding_test.cpp
)

target_compile_definitions(
  MODULE_CONTEXT_binding_test
  PRIVATE
  ABACUS_CAN_BIND_SIMULATION_CONTEXT
  ABACUS_CAN_READ_SIMULATION_CONTEXT)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(
  NAME MODULE_CONTEXT_mapping_manifest_test
  COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/check_context_mapping.py
          --repo-root ${PROJECT_SOURCE_DIR})
