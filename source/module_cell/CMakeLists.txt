add_subdirectory(module_symmetry)
add_subdirectory(module_neighbor)
if(ENABLE_PAW)
  add_subdirectory(module_paw)
endif()

add_library(
    cell
    OBJECT
    atom_pseudo.cpp
    atom_spec.cpp
    pseudo_nc.cpp
    read_pp.cpp
    read_pp_upf100.cpp
    read_pp_upf201.cpp
    read_pp_blps.cpp
    read_pp_vwr.cpp
    unitcell.cpp
    read_atoms.cpp
    read_cell_pseudopots.cpp
    setup_nonlocal.cpp
    klist.cpp
    parallel_kpoints.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(cell)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
  add_subdirectory(test_pw)
endif()
