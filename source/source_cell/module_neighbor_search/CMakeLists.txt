add_library(
    neighbor_search
    OBJECT
    bin_manager.cpp
    neighbor_search.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(neighbor_search)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()