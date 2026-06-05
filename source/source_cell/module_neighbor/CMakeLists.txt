add_library(
    neighbor
    OBJECT
    sltk_atom.cpp
    sltk_atom_arrange.cpp
    sltk_grid.cpp
    sltk_grid_driver.cpp
    atom_pack.cpp
    mpi_domain.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(neighbor)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
