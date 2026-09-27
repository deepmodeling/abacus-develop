add_library(librpa_io INTERFACE)
target_include_directories(librpa_io INTERFACE
  ${ABACUS_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(test)
endif()
