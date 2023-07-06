remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET cell_neighbor_sltk_adjacent_set
  SOURCES sltk_adjacent_set_test.cpp ../sltk_adjacent_set.cpp
)

AddTest(
  TARGET cell_neighbor_sltk_atom
  SOURCES sltk_atom_test.cpp ../sltk_atom.cpp ../sltk_adjacent_set.cpp
)

AddTest(
  TARGET cell_neighbor_sltk_atom_arrange
  LIBS MPI::MPI_CXX
  SOURCES ./sltk_atom_arrange_test.cpp  ../sltk_atom_arrange.cpp 
)


AddTest(
  TARGET cell_neighbor_sltk_grid
  LIBS MPI::MPI_CXX
  SOURCES ./sltk_grid_test.cpp 
)