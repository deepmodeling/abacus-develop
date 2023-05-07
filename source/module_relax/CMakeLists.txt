add_library(
    relax
    OBJECT

    relax_driver.cpp

    relax_new/relax.cpp
    relax_new/line_search.cpp

    relax_old/relax_old.cpp
    relax_old/bfgs_basic.cpp
    relax_old/ions_move_basic.cpp
    relax_old/ions_move_bfgs.cpp
    relax_old/ions_move_cg.cpp
    relax_old/ions_move_sd.cpp
    relax_old/ions_move_methods.cpp
    relax_old/lattice_change_basic.cpp
    relax_old/lattice_change_cg.cpp
    relax_old/lattice_change_methods.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(relax)
endif()

if(BUILD_TESTING)
  add_subdirectory(relax_new/test)
endif()

if(BUILD_TESTING)
  add_subdirectory(relax_old/test)
endif()
