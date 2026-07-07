set(neighbor_search_sources
    bin_manager.cpp
    domain_decomposition.cpp
    neighbor_search.cpp
    page_allocator.cpp
    unitcell_lite.cpp
)

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
