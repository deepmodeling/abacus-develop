set(neighbor_search_sources
    bin_manager.cpp
    neighbor_search.cpp
    page_allocator.cpp
)

if(USE_CUDA AND ENABLE_GPU_NEIGHBOR_LIST)
  list(APPEND neighbor_search_sources gpu_neighbor_list.cpp)
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
