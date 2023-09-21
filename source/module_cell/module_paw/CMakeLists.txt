add_library(
    paw
    OBJECT
    paw_element.cpp
    paw_sphbes.cpp
    paw_cell.cpp
    paw_cell_libpaw.cpp
    paw_atom.cpp
)

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
 target_link_libraries(paw libpaw_interface -lifcore)
else()
 target_link_libraries(paw libpaw_interface)
endif()
add_dependencies(paw libpaw_interface)

if(ENABLE_COVERAGE)
  add_coverage(paw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()