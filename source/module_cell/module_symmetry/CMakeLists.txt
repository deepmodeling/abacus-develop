add_library(
    symmetry
    OBJECT
    symm_other.cpp
    symmetry_basic.cpp
    symmetry.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(symmetry)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
