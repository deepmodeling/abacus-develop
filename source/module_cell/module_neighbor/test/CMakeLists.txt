remove_definitions(-D__DEEPKS)
remove_definitions(-D__CUDA)
remove_definitions(-D__ROCM)
remove_definitions(-D__EXX)

AddTest(
  TARGET cell_neighbor_sltk_adjacent_set
  SOURCES ../sltk_adjacent_set.cpp ./sltk_adjacent_set_test.cpp
)

AddTest(
  TARGET cell_neighbor_load_atom
  SOURCES ../sltk_atom_input.cpp ./sltk_atom_input_test.cpp
)

