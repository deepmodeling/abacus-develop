add_library(
    paw
    OBJECT
    paw_element.cpp
    paw_sphbes.cpp
    paw_cell.cpp
    paw_cell_libpaw.cpp
    paw_atom.cpp
)

if(ENABLE_PAW)
  target_link_libraries(paw libpaw_interface -lifcore)
  add_dependencies(paw libpaw_interface)
endif()

if(ENABLE_COVERAGE)
  add_coverage(paw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()