add_library(
    relax
    OBJECT

    relax_driver.cpp

    relax_sync.cpp
    line_search.cpp
    
    bfgs.cpp
    lbfgs.cpp
    relax_nsync.cpp
    bfgs_basic.cpp
    ions_move_basic.cpp
    ions_move_bfgs.cpp
    ions_move_cg.cpp
    ions_move_sd.cpp
    ions_move_methods.cpp
    lattice_change_basic.cpp
    lattice_change_cg.cpp
    lattice_change_methods.cpp
    matrix_methods.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(relax)
endif()

if(BUILD_TESTING)
  if(ENABLE_MPI)
    add_subdirectory(test)
  endif()
   
endif()
