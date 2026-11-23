add_library(librpa_io OBJECT librpa_bz_sampling.cpp)
target_include_directories(librpa_io PUBLIC ${ABACUS_SOURCE_DIR}
                                            ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_TESTING AND ENABLE_MPI)
  add_subdirectory(unittests)
endif()
