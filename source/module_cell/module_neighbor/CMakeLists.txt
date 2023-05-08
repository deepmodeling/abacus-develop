add_library(
    neighbor
    OBJECT
    sltk_adjacent_set.cpp
    sltk_atom.cpp
    sltk_atom_arrange.cpp
    sltk_atom_input.cpp
    sltk_grid.cpp
    sltk_grid_driver.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(neighbor)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()