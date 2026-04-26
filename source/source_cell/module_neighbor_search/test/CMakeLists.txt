remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)



AddTest(
  TARGET MODULE_CELL_NEIGHBOR_neighbor_search
  LIBS parameter ${math_libs} base device 
  SOURCES
    test_neighbor_search.cpp
    ../neighbor_search.cpp
    ../bin_manager.cpp
)


