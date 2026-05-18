add_library(
    neighbor
    OBJECT
    sltk_atom.cpp
    sltk_atom_arrange.cpp
    sltk_grid.cpp
    sltk_grid_driver.cpp
    neighlist_adapter.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(neighbor)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
endif()
