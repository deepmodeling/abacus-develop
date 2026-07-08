abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)
abacus_disable_feature_definitions(__EXX)



# Split tests: create one test target per test source to keep failures isolated.

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_neighbor_search
  LIBS parameter base device
  SOURCES
    neighbor_search_test.cpp
  ../neighbor_search.cpp
  ../bin_manager.cpp
  ../page_allocator.cpp
  ../unitcell_lite.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_bin_manager
  LIBS parameter base device
  SOURCES
    bin_manager_test.cpp
    ../bin_manager.cpp
    ../page_allocator.cpp
)

AddTest(
  TARGET MODULE_CELL_NEIGHBOR_allocator_and_list
  LIBS parameter base device
  SOURCES
  neighbor_list_test.cpp
  ../page_allocator.cpp
)


