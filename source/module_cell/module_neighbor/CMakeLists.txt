add_library(
    neighbor
    OBJECT
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
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()