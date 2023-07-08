add_library(
    paw
    OBJECT
    paw_element.cpp
    paw_interface.cpp
    paw_sphbes.cpp
    paw_cell.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(paw)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()