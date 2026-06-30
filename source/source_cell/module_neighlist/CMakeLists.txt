set(neighbor_search_sources
    bin_manager.cpp
    neighbor_search.cpp
    page_allocator.cpp
    unitcell_lite.cpp
)

if(ENABLE_MPI)
    list(APPEND neighbor_search_sources domain_decomposition.cpp)
endif()

add_library(
    neighbor_search
    OBJECT
    ${neighbor_search_sources}
)

if(ENABLE_COVERAGE)
  add_coverage(neighbor_search)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
