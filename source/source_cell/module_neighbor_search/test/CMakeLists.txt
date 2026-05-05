remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)



# Split tests: create one test target per test source to keep failures isolated.

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_neighbor_search
  LIBS parameter ${math_libs} base device
  SOURCES
    neighbor_search_test.cpp
  ../neighbor_search.cpp
  ../bin_manager.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_bin_manager
  LIBS parameter ${math_libs} base device
  SOURCES
    bin_manager_test.cpp
    ../bin_manager.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_allocator_and_list
  LIBS parameter ${math_libs} base device
  SOURCES
  neighbor_list_test.cpp
)


